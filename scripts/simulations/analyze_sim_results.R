library(tidyverse)
library(arrow)
library(glue)
library(patchwork)

source('scripts/helpers.R') 

### Analysis of Standardization Metric Summary Statistics
df_sm <- 
  bind_rows(
    map_dfr(dir('sim_outputs/n100/standardization_metrics/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n500/standardization_metrics/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n1000/standardization_metrics/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n5000/standardization_metrics/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n10000/standardization_metrics/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n25000/standardization_metrics/', full.names = T), read_parquet)
  )

true_sigma <- read_csv('sim_inputs/true_sigma_ratios.csv') 


### plot Mean Sigma 2 Ratio
df_summary <- 
  df_sm %>% 
  group_by(shift, trt_change, n_subjects, estimator, trial_id) %>% 
  summarise('mean_sigma2_ratio' = mean(sigma2_ratio),
            'sd_sigma2_ratio' = sd(sigma2_ratio))

df_cross <- 
  df_summary %>% 
  filter(estimator == 'chi_cross') %>% 
  mutate('n_subjects' = paste('n =', gsub('000$', ',000', n_subjects))) %>% 
  bind_rows(true_sigma %>% 
              filter(estimator == 'chi_DRj') %>% 
              mutate('n_subjects' = 'Truth') %>% 
              rename('mean_sigma2_ratio' = sigma2_ratio)) %>% 
  mutate('shift' = case_when(shift == 'No Covariate Shift' ~ "No\nCovariate Shift",
                             shift == 'Linear Covariate Shift' ~ "Linear\nCovariate Shift",
                             shift == 'Flexible Covariate Shift' ~ "Flexible\nCovariate Shift")) %>% 
  mutate('shift' = fct_rev(shift)) %>% 
  mutate('n_subjects' = as.factor(paste0('"', n_subjects, '"'))) %>% 
  mutate('n_subjects' = fct_relevel(as.factor(n_subjects), 
                                    paste0('"', c('n = 100', 'n = 500', 'n = 1,000', 'n = 5,000', 
                                                  'n = 10,000', 'n = 25,000', 'Truth'), '"')))

p1 <- 
  ggplot(df_cross, aes(x = trial_id, y = mean_sigma2_ratio)) + 
  facet_grid(fct_rev(paste0('"', shift, '"')) ~ n_subjects, labeller = label_parsed) +
  geom_point(aes(color = trt_change), size = 3, alpha = 0.5) + 
  scale_y_continuous(labels = scales::percent) + 
  theme(legend.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 16)) + 
  labs(x = 'Trial Index (m)',
       y = expression(paste('Average ', hat(theta)[m])),
       color = 'Underlying Treatment\nEffect Structure',
       title = expression(paste('Estimated Proportion of Variation in ', widehat(chi)['j,m'],  ' Not Attributable to Distribution Shift')),
       tag = 'A)'
  )


### Analysis of Loss Fx
df_loss <- 
  bind_rows(
    map_dfr(dir('sim_outputs/n100/loss_fx/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n500/loss_fx/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n1000/loss_fx/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n5000/loss_fx/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n10000/loss_fx/', full.names = T), read_parquet),
    map_dfr(dir('sim_outputs/n25000/loss_fx/', full.names = T), read_parquet),
  )

### Use SD-based difference Rule
tolerance_sd <- function(c) {
  df_tolerance_sd <- 
    df_loss %>% 
    mutate('complexity' = case_when(projection == 'constant' ~ 1,
                                    projection == 'linear' ~ 2,
                                    projection == 'cubic' ~ 3,
                                    projection == 'spline_2' ~ 4,
                                    projection == 'spline_3' ~ 5)) %>% 
    group_by(sim_id, trt_change, shift, n_subjects) %>% 
    mutate('min_risk_adj' = min(risk_adj)) %>% 
    mutate('tolerance_sd' =  weighted_sd[which.min(risk_adj)]) %>% 
    mutate('within_tolerance_sd' = abs(risk_adj - min_risk_adj)/tolerance_sd <= c) %>% 
    mutate('best_within_tol_sd' = complexity == min(complexity[within_tolerance_sd])) %>% 
    group_by(parameter_id, trt_change, shift, n_subjects, projection) %>% 
    summarise('best_proj_within_tol_sd' = mean(best_within_tol_sd)) %>% 
    ungroup() %>% 
    mutate('subjects_f' = fct_reorder(as.factor(prettyNum(n_subjects, big.mark = ',')), n_subjects)) %>% 
    mutate('shift' = case_when(shift == 'No Covariate Shift' ~ "No\nCovariate Shift",
                               shift == 'Linear Covariate Shift' ~ "Linear\nCovariate Shift",
                               shift == 'Flexible Covariate Shift' ~ "Flexible\nCovariate Shift")) %>% 
    mutate('shift' = fct_rev(shift)) %>% 
    mutate('projection' = gsub('psi_', '', projection),
           'projection' = gsub('_', ' (df = ', tools::toTitleCase(projection)), 
           'projection' = ifelse(grepl('Spline', projection), paste0(projection, ')'), projection),
           'projection' = case_when(projection == 'Spline (df = 2)' ~ 'Spline (2 Knots)',
                                    projection == 'Spline (df = 3)' ~ 'Spline (3 Knots)',
                                    T ~ projection),
           'projection' = fct_relevel(projection, 'Constant', 'Linear', 'Cubic', 'Spline (2 Knots)', 'Spline (3 Knots)')) %>% 
    mutate('c' = c)
  
}


df_tolerance_sd <- 
  map_dfr(c(0, 0.25, 0.5, 0.75, 1), tolerance_sd) %>% 
  mutate('trt_change' = gsub('\\+\\s', '\\+\n', trt_change)) %>% 
  mutate('trt_change' = fct_relevel(paste0('"', trt_change, '"'), paste0('"', levels(as.factor(trt_change)), '"')))


p2 <- 
  ggplot(df_tolerance_sd, aes(x = subjects_f, y = best_proj_within_tol_sd)) +
  facet_grid(fct_rev(paste0('"', shift, '"')) ~ trt_change,
             labeller = label_parsed) +
  geom_point(aes(color = projection, shape = as.factor(c)), size = 5, alpha = 0.5) + 
  scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
  scale_shape_manual(values = c(1,16,6,3,8))  +
  theme(legend.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16),
        legend.box = 'vertical') + 
  labs(x = '# of Subjects',
       y = expression(paste('% of Simulations Minimizing Risk (Within ', c, epsilon, ')')),
       color = expression(paste('Projection ', psi(m * ";" * bold(beta)))),
       shape = '# of Standard Deviations Tolerance (c)',
       title = 'Frequency of Projection Selection',
       tag = 'B)'
       # subtitle = 'Simulation Results'
  )

p2/p1
ggsave('figures/sim_results.pdf', height = 2 * 7, width = 12)

### Summary of scenarios 
input_summary <- 
  read_csv('sim_inputs/input_summary.csv') %>% 
  mutate('parameter_id' = 1:nrow(.),
         'region' = case_when(outcome_model %in% c('constant', 'effect_mod') ~ 'Null 0',
                              outcome_model %in% c('spline', 'linear') ~ 'Null 1',
                              covariate_shift == 'no_covariate_shift' & outcome_model %in% c('linear_effect_mod', 'spline_effect_mod') ~ 'Null 1',
                              T ~ 'Alternative'))

df_hypothesis <- 
  df_hypothesis %>% 
  inner_join(input_summary, by = 'parameter_id')

### Plot Correct Rate 
compute_delta_stats <- function(delta) {
  df_summary <- 
    df_hypothesis %>% 
    group_by(shift, trt_change, region, n_subjects) %>% 
    summarise('delta' = delta, 
              'reject0' = mean(qlow > delta),
              'reject1' = mean(qhigh < 1-delta),
              'reject_both' = mean(qlow > delta & qhigh < 1-delta)) %>% 
    mutate('correct_rate' = 
             case_when(region == 'Null 0' ~ 1 - reject0,
                       region == 'Null 1' ~ 1 - reject1,
                       region == 'Alternative' ~ reject_both),
           'error_rate' = 
             case_when(region != 'Alternative' ~ reject_both,
                       region == 'Alternative' ~ 1-reject_both)
    ) %>% 
    ungroup() %>% 
    mutate('shift' = case_when(shift == 'No Covariate Shift' ~ "No\nCovariate Shift",
                               shift == 'Linear Covariate Shift' ~ "Linear\nCovariate Shift",
                               shift == 'Flexible Covariate Shift' ~ "Flexible\nCovariate Shift")) 
  
  
  return(df_summary)
}

df_delta <- 
  map_dfr(c(0.01, 0.025, 0.05, 0.075, 0.1), compute_delta_stats) %>% 
  mutate('subjects_f' = fct_reorder(as.factor(prettyNum(n_subjects, big.mark = ',')), n_subjects)) %>% 
  mutate('trt_change' = gsub('\\+\\s', '\\+\n', trt_change)) %>% 
  mutate('trt_change' = fct_relevel(paste0('"', trt_change, '"'), paste0('"', levels(as.factor(trt_change)), '"')))


ggplot(df_delta, aes(x = subjects_f, y = correct_rate)) + 
  facet_grid(fct_rev(paste0('"', shift, '"')) ~ trt_change,
             labeller = label_parsed) +
  geom_point(aes(shape = region, color = as.factor(delta)), size = 3, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent) + 
  scale_shape_manual(
    breaks = c('Null 0', 'Null 1', 'Alternative'), 
    values = c(16, 17, 8),
    labels = list(
      expression(paste('[0,', delta, ']')),
      expression(paste('[1-' , delta,  ', 1]')),
      expression(paste('(', delta, ', 1-', delta, ')'))
    )
  ) +
  theme(legend.title = element_text(size = 15),
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box = 'vertical') + 
  labs(x = '# of Subjects',
       y = '% of Simulations Placing Variability Ratio in\nCorrect Hypothesis Region',
       color = expression(paste('Boundary Margin (', delta, ')')),
       title = expression(paste('Hypothesis Testing for ', theta)),
       # subtitle = 'Simulation Results',
       shape = 'True Hypothesis Region')
ggsave('figures/hypothesis_testing_results.png', width = 16, height = 9)
ggsave('figures/hypothesis_testing_results.pdf', width = 12, height = 7)
