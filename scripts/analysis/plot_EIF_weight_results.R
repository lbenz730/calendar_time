library(tidyverse)
library(glue)
library(patchwork)

source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Results
df_projections <- read_csv(glue('{data_dir}/tv_effects/weight_projections.csv'))

df_plot <- 
  df_projections %>% 
  mutate('outcome' = case_when(outcome == 'delta_6mo' ~ "'6 Months'",
                               outcome == 'delta_1yr' ~ "'1 Year'",
                               outcome == 'delta_2yr' ~ "'2 Years'", 
                               outcome == 'delta_3yr' ~ "'3 Years'"),
         'outcome' = fct_relevel(outcome, "'6 Months'"), 
         'comparison' = case_when(procedures == 'AGB/RYGB/SLEEVE' ~ "'Surgery vs.\nNo Surgery'",
                                  procedures == 'RYGB' ~ "'RYGB vs.\nNo Surgery'",
                                  procedures == 'SLEEVE' ~ "'SG vs.\nNo Surgery'",
                                  procedures == 'RYGB/SLEEVE' ~ "'RYGB vs.\nSG'"),
         'comparison' = fct_relevel(comparison,
                                    "'Surgery vs.\nNo Surgery'", 
                                    "'RYGB vs.\nNo Surgery'", 
                                    "'SG vs.\nNo Surgery'",
                                    "'RYGB vs.\nSG'"),
         'projection' = gsub('_', ' (df = ', tools::toTitleCase(projection)), 
         'projection' = ifelse(grepl('Spline', projection), paste0(projection, ')'), projection)) %>% 
  mutate('projection' = case_when(projection == 'Spline (df = 2)' ~ 'Spline (2 Knots)',
                                  projection == 'Spline (df = 3)' ~ 'Spline (3 Knots)',
                                  T ~ projection)) %>% 
  mutate('projection' = fct_relevel(projection, 'Constant', 'Linear', 'Cubic', 'Spline (2 Knots)', 'Spline (3 Knots)'))

p1 <- 
  ggplot(df_plot, aes(x = trial_id, y = chi_DR)) + 
  facet_grid(comparison ~ outcome, scales = 'free_y', labeller = label_parsed) +  
  geom_point(data = filter(df_plot, weights == 'uniform', projection == 'Constant'),  alpha = 0.2, size = 3) + 
  geom_ribbon(aes(ymin = psi_hat + qnorm(0.025) * sd,
                  ymax = psi_hat + qnorm(0.975) * sd,
                  fill = projection), alpha = 0.2) +
  geom_line(aes(y = psi_hat, color = projection), lwd = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Calendar Time-Specific Treatment Effect Estimates',
       x = 'Trial Index (m)',
       y = 'Difference in Mean %-Weight Change',
       plot.tag = 'A)',
       color = expression(psi(m * ";" * bold(beta))),
       fill = expression(psi(m * ";" * bold(beta)))) +
  theme(plot.tag = element_text(size = 16))

df_transport <- 
  map_dfr(1:84, ~{
    read_csv(glue('{data_dir}/tv_effects/weight_DR_results_baseline_trial_{.x}.csv')) %>% 
      mutate('baseline_trial' = .x)
  }) %>% 
  mutate('outcome' = case_when(outcome == 'delta_6mo' ~ '6 Months',
                               outcome == 'delta_1yr' ~ '1 Year',
                               outcome == 'delta_2yr' ~ '2 Years', 
                               outcome == 'delta_3yr' ~ '3 Years'),
         'outcome' = fct_relevel(outcome, '6 Months'), 
         'comparison' = case_when(procedures == 'AGB/RYGB/SLEEVE' ~ 'Surgery vs. No Surgery',
                                  procedures == 'RYGB' ~ 'RYGB vs. No Surgery',
                                  procedures == 'SLEEVE' ~ 'SG vs. No Surgery',
                                  procedures == 'RYGB/SLEEVE' ~ 'RYGB vs. SG'),
         'comparison' = fct_relevel(comparison,
                                    'Surgery vs. No Surgery', 
                                    'RYGB vs. No Surgery', 
                                    'SG vs. No Surgery',
                                    'RYGB vs. SG'))

p2 <- 
  ggplot(df_transport %>% filter(comparison %in% c('Surgery vs. No Surgery', 
                                                   'SG vs. No Surgery')), aes(x = trial_id, y = chi_cross)) + 
  facet_grid(comparison~outcome, scales = 'free_y') +  
  geom_point(aes(color = baseline_trial)) + 
  scale_y_continuous(labels = scales::percent)  +
  scale_color_viridis_c(option = 'H') + 
  theme(plot.tag = element_text(size = 16)) +
  labs(x = 'Treatment Initaition Time Trial Index (m)',
       y = 'Difference in Mean %-Weight Change',
       title = expression(paste('Illustration of ', widehat(chi)["j,m"], ' for Select Comparisons')),
       color = 'Covariate Distribution Trial Index (j)',
       plot.tag = 'B)')

p1/p2

ggsave('figures/final_analysis_figure.pdf', height = 2 * 7, width = 12)


### Hypothesis Testing for Sigma2 Ratio
df_hypothesis <- 
  map_dfr(1:16, ~read_csv(glue('{data_dir}/tv_effects/EIF_weight_hypothesis_tests_{.x}.csv'))) 

df_hypothesis <- 
  df_hypothesis %>% 
  mutate('outcome' = case_when(outcome == 'delta_6mo' ~ '6 Months',
                               outcome == 'delta_1yr' ~ '1 Year',
                               outcome == 'delta_2yr' ~ '2 Years', 
                               outcome == 'delta_3yr' ~ '3 Years'),
         'outcome' = fct_relevel(outcome, '6 Months'), 
         'comparison' = case_when(procedures == 'AGB/RYGB/SLEEVE' ~ 'Surgery vs. No Surgery',
                                  procedures == 'RYGB' ~ 'RYGB vs. No Surgery',
                                  procedures == 'SLEEVE' ~ 'SG vs. No Surgery',
                                  procedures == 'RYGB/SLEEVE' ~ 'RYGB vs. SG'),
         'comparison' = fct_relevel(comparison,
                                    'Surgery vs. No Surgery', 
                                    'RYGB vs. No Surgery', 
                                    'SG vs. No Surgery',
                                    'RYGB vs. SG'),
         'interval' = paste0('(', sprintf('%0.3f', qlow), ', ', sprintf('%0.3f', qhigh), ')'))


df_risk <- read_csv(glue('{data_dir}/tv_effects/weight_projection_risk.csv'))
df_risk <- 
  df_risk %>% 
  mutate('outcome' = case_when(outcome == 'delta_6mo' ~ '6 Months',
                               outcome == 'delta_1yr' ~ '1 Year',
                               outcome == 'delta_2yr' ~ '2 Years', 
                               outcome == 'delta_3yr' ~ '3 Years'),
         'outcome' = fct_relevel(outcome, '6 Months'), 
         'comparison' = case_when(procedures == 'AGB/RYGB/SLEEVE' ~ 'Surgery vs. No Surgery',
                                  procedures == 'RYGB' ~ 'RYGB vs. No Surgery',
                                  procedures == 'SLEEVE' ~ 'SG vs. No Surgery',
                                  procedures == 'RYGB/SLEEVE' ~ 'RYGB vs. SG'),
         'comparison' = fct_relevel(comparison,
                                    'Surgery vs. No Surgery', 
                                    'RYGB vs. No Surgery', 
                                    'SG vs. No Surgery',
                                    'RYGB vs. SG')) %>% 
  mutate('projection' =  gsub('_', ' (df = ', tools::toTitleCase(gsub('psi_', '', projection))), 
         'projection' = ifelse(grepl('Spline', projection), paste0(projection, ')'), projection)) %>% 
  mutate('projection' = fct_relevel(projection, 'Constant', 'Linear', 'Cubic', 'Spline (df = 2)', 'Spline (df = 3)')) %>% 
  mutate('projection' = case_when(projection == 'Spline (df = 2)' ~ 'Spline (2 Knots)',
                                  projection == 'Spline (df = 3)' ~ 'Spline (3 Knots)',
                                  T ~ projection)) %>% 
  mutate('projection' = fct_relevel(projection, 'Constant', 'Linear', 'Cubic', 'Spline (2 Knots)', 'Spline (3 Knots)')) %>% 
  group_by(outcome, procedures, weights) %>% 
  mutate('min_risk' = risk == min(risk),
         'min_risk_adj' = risk_adj == min(risk_adj),
         'within_tol' = abs(risk_adj-min(risk_adj))/weighted_sd[which.min(risk_adj)] <= 1/4,
         'min_risk_within_tol' = as.numeric(projection) == min(as.numeric(projection[within_tol]))) %>% 
  select(-within_tol) %>% 
  ungroup()

### Table S1
df_risk %>% 
  group_by(comparison, outcome) %>% 
  mutate('sd_diff' = abs( risk_adj - min(risk_adj) )/weighted_sd[which.min(risk_adj)]) %>% 
  mutate('summary' = paste0(sprintf('%0.3f', risk_adj), ' (', sprintf('%0.2f', sd_diff), ')')) %>% 
  ungroup() %>% 
  select(comparison, outcome, projection, summary) %>% 
  pivot_wider(names_from = projection,
              values_from = summary) %>% 
  arrange(comparison, outcome) %>% 
  kbl(align = 'c', format = 'latex') %>% 
  collapse_rows(1)



### Latex Table for Main Paper
df_latex <- 
  df_hypothesis %>% 
  select(comparison, outcome, interval) %>% 
  inner_join( 
    df_risk %>% 
      group_by(comparison, outcome) %>% 
      summarise('minimizer' = projection[min_risk_adj],
                'tol_minimizer' = projection[min_risk_within_tol]),
    by = c('comparison', 'outcome')) %>% 
  select(comparison, outcome, minimizer, tol_minimizer, interval) %>% 
  arrange(comparison, outcome) 

kbl(df_latex, align = 'c', format = 'latex') %>% 
  collapse_rows(1)

