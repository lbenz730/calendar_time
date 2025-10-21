library(tidyverse)
library(glue)
library(arrow)
library(patchwork)
source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R')

set.seed(04071318)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 

summary <- 
  df_trials %>% 
  filter(eligible) %>% 
  group_by(trial_id) %>% 
  summarise('n_control' = sum(1-surgery),
            'n_surg' = sum(surgery),
            'n_rygb' = sum(bs_type == 'RYGB'),
            'n_sg' = sum(bs_type == 'SLEEVE')) %>% 
  pivot_longer(cols = contains('n_'),
               names_to = 'group',
               values_to = 'n',
               names_prefix = 'n_') %>% 
  mutate('group' = case_when(group == 'rygb' ~ 'RYGB',
                             group == 'sg' ~ 'SG',
                             group == 'surg' ~ 'Surgery',
                             group == 'control' ~ 'No Surgery'))

ggplot(summary, aes(x = trial_id, y = n)) + 
  facet_wrap(~ifelse(group == 'No Surgery', 'No Surgery', 'Surgery'), scales = 'free_y') + 
  geom_point(aes(color = group)) + 
  geom_line(aes(color = group)) + 
  scale_y_continuous(labels = ~scales::number(.x, big.mark = ',')) + 
  labs(x = 'Trial Index (m)',
       y = '# of Eligible Patients',
       color = '')

ggsave('figures/population_size.png', height = 9/1.2, width = 16/1.2)


df_trials <- 
  df_trials %>% 
  filter(eligible) %>% 
  filter(!is.na(delta_3yr))

### Split datasets into subgroups
df_control <- 
  df_trials %>% 
  filter(surgery == 0)

df_surgery <- 
  df_trials %>% 
  filter(surgery == 1) 

df_rygb <- 
  df_surgery %>% 
  filter(bs_type == 'RYGB') 

df_sleeve <- 
  df_surgery %>% 
  filter(bs_type == 'SLEEVE') 




plot_dists <- function(df, title, subtitle) {
  ### Baseline BMI
  df_bmi <- 
    df %>% 
    group_by(trial_id) %>% 
    summarise('mean_bmi' = mean(baseline_bmi))
  
  plot_bmi <- 
    ggplot(df_bmi, aes(x = trial_id, y = mean_bmi)) + 
    facet_wrap(~'BMI') +
    geom_line() +
    geom_point(size = 2) +
    geom_smooth() + 
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    labs(x = 'Trial Index (m)',
         y = '')
  
  ### Baseline Age
  df_age <- 
    df %>% 
    group_by(trial_id) %>% 
    summarise('mean_age' = mean(baseline_age))
  
  plot_age <- 
    ggplot(df_age, aes(x = trial_id, y = mean_age)) + 
    facet_wrap(~'Age') +
    geom_line() +
    geom_point(size = 2) +
    geom_smooth() + 
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    labs(x = 'Trial Index (m)',
         y = '')
  
  ### Sex
  df_sex <- 
    df %>% 
    group_by(trial_id, gender = as.factor(gender), .drop = F) %>% 
    count() %>% 
    group_by(trial_id) %>% 
    mutate('pct' = n/sum(n))
  
  plot_sex <-
    ggplot(df_sex, aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Sex') +
    geom_line(aes(color = gender)) +
    geom_point(aes(color = gender), size = 2) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    scale_color_manual(values = c('hotpink', 'dodgerblue')) + 
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  ### Site
  df_site <- 
    df %>% 
    group_by(trial_id, site = as.factor(site), .drop = F) %>% 
    count() %>% 
    group_by(trial_id) %>% 
    mutate('pct' = n/sum(n),
           'site' = as.character(site)) %>% 
    mutate('site' = case_when(site == 'GH' ~ 'KPWA', 
                              T ~ paste0('KP', site)))
  
  plot_site <-
    ggplot(df_site, aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Site') +
    geom_line(aes(color = site)) +
    geom_point(aes(color = site)) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    scale_color_manual(values = c('darkgreen', 'red', 'navy')) + 
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  ### Race
  df_race <- 
    df %>% 
    group_by(trial_id, race = as.factor(race), .drop = F) %>% 
    count() %>% 
    group_by(trial_id) %>% 
    mutate('pct' = n/sum(n),
           'race' = as.character(race)) %>% 
    mutate('race' = case_when(race == 'BA' ~ 'Black', 
                              race == 'WH' ~ 'White', 
                              T ~ 'Other')) %>% 
    mutate('race' = fct_relevel(race, 'Black', 'White', 'Other'))
  
  plot_race <-
    ggplot(df_race, aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Race') +
    geom_line(aes(color = race)) +
    geom_point(aes(color = race), size = 2) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  ### Hypertension + Dyslipidemia + Medications
  df_ht <- 
    df %>% 
    group_by(trial_id) %>% 
    summarise('Hypertension' = mean(hypertension),
              'Hypertension Rx' = mean(hypertension_rx)) %>% 
    pivot_longer(cols = -trial_id,
                 names_to = 'covariate',
                 values_to = 'pct')
  
  df_dys <- 
    df %>% 
    group_by(trial_id) %>% 
    summarise('Dyslipidemia' = mean(dyslipidemia),
              'Antilipemic Rx' = mean(antilipemic_rx)) %>% 
    pivot_longer(cols = -trial_id,
                 names_to = 'covariate',
                 values_to = 'pct')
  
  plot_dys <- 
    bind_rows(df_dys, df_ht) %>% 
    mutate('covariate' = fct_relevel(covariate, 'Dyslipidemia', 'Antilipemic Rx')) %>% 
    ggplot(aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Hypertension + Dyslipidemia') +
    geom_line(aes(color = covariate)) +
    geom_point(size = 2, aes(color = covariate)) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    scale_color_manual(values = c('navy', 'skyblue', 'red', 'light coral')) +
    guides(color = guide_legend(nrow = 2)) +
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  ### T2DM 
  df_diabetes <- 
    df %>% 
    group_by(trial_id) %>% 
    summarise('T2DM' = mean(t2dm),
              'Insulin Rx' = mean(insulin)) %>% 
    pivot_longer(cols = -trial_id,
                 names_to = 'covariate',
                 values_to = 'pct') %>% 
    mutate('covariate' = fct_relevel(covariate, 'T2DM')) 
  
  plot_t2dm <- 
    ggplot(df_diabetes, aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Diabetes') +
    geom_line(aes(color = covariate)) +
    geom_point(aes(color = covariate), size = 2) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    scale_color_manual(values = c('purple', 'plum')) +
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  df_smoke <- 
    df %>% 
    group_by(trial_id, smoking_status = as.factor(smoking_status), .drop = F) %>% 
    count() %>% 
    group_by(trial_id) %>% 
    mutate('pct' = n/sum(n)) %>% 
    mutate('smoking_status' = as.character(smoking_status)) %>% 
    mutate('smoking_status' = case_when(smoking_status == 'no_self_report' ~ 'N/A',
                                        T ~ tools::toTitleCase(smoking_status))) %>% 
    mutate('smoking_status' = fct_relevel(smoking_status, 'Current', 'Former', 'Never', 'N/A'))
  
  plot_smoke <-
    ggplot(df_smoke, aes(x = trial_id, y = pct)) + 
    facet_wrap(~'Self-Reported Smoking Status') +
    geom_line(aes(color = smoking_status)) +
    geom_point(aes(color = smoking_status), size = 2) +
    scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) + 
    scale_color_manual(values = c('black', 'brown', 'grey', 'skyblue')) +
    labs(x = 'Trial Index (m)', 
         y = '',
         color = '')
  
  
  final_plot <- 
    (plot_age + plot_bmi + plot_dys + plot_t2dm + plot_race + plot_sex + plot_smoke + plot_site) + 
    plot_layout(ncol = 4) + 
    plot_layout(axis_titles = 'collect') + 
    plot_annotation(title = title,
                    subtitle = subtitle)
  
  
  
  return(final_plot)
  
}


### Make Plots
plot_dists(df = df_surgery,
           title = 'Distribution of Key Covariates in Study Population',
           subtitle = 'Eligible Patients Undergoing Bariatric Surgery')
ggsave('figures/confounder_dist/confounder_dist_surgery.png', height = 9, width = 16)
ggsave('figures/confounder_dist/confounder_dist_surgery.pdf', height = 9, width = 16)

plot_dists(df = df_control,
           title = 'Distribution of Key Covariates in Study Population',
           subtitle = 'Eligible Patients Not Undergoing Bariatric Surgery')
ggsave('figures/confounder_dist/confounder_dist_control.png', height = 9, width = 16)

plot_dists(df = df_rygb,
           title = 'Distribution of Key Covariates in Study Population',
           subtitle = 'Eligible Patients Undergoing Roux-en-Y Gastric Bypass')
ggsave('figures/confounder_dist_rygb.png', height = 9, width = 16)

plot_dists(df = df_sleeve,
           title = 'Distribution of Key Covariates in Study Population',
           subtitle = 'Eligible Patients Undergoing Sleeve Gastrectomy')
ggsave('figures/confounder_dist/confounder_dist_sg.png', height = 9, width = 16)

