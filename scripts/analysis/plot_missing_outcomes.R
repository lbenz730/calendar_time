library(tidyverse)
library(glue)
library(arrow)
source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 


### Summary of Missing Data Over Time
df_na <- 
  df_trials %>% 
  filter(eligible == T) %>% 
  group_by(trial_id, surgery) %>% 
  summarise('pct_06mo' = mean(!is.na(delta_6mo)),
            'pct_1yr' = mean(!is.na(delta_1yr)),
            'pct_2yr' = mean(!is.na(delta_2yr)),
            'pct_3yr' = mean(!is.na(delta_3yr))) %>% 
  ungroup() %>% 
  pivot_longer(cols = contains('pct_'),
               names_to = 'time_frame',
               values_to = 'pct_na',
               names_prefix = 'pct_') %>% 
  mutate('time_frame' = case_when(time_frame == '06mo' ~ '6 Months',
                                  time_frame == '1yr' ~ '1 Year',
                                  time_frame == '2yr' ~ '2 Years',
                                  time_frame == '3yr' ~ '3 Years')) %>% 
  mutate('time_frame' = fct_relevel(time_frame, '6 Months')) %>% 
  mutate('surgery' = ifelse(surgery == 0, 'No Surgery', 'Surgery')) 


ggplot(df_na, aes(x = trial_id, y = pct_na)) + 
  facet_wrap(~time_frame) + 
  geom_point(aes(color = as.factor(surgery))) +
  geom_line(aes(color = as.factor(surgery))) + 
  scale_x_continuous(limits = c(0, 84), breaks = seq(0, 84, 12)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x = 'Trial Index (m)',
       y = '% of Eligible Patients with Observed Outcome',
       color = '',
       title = 'Outcome Ascertainment Rate by Surgical Status Over Time',
       subtitle = 'Eligible Patients in DURABLE Database (2005-2011)')
ggsave('figures/sensitiviy/outcome_missing_data.png', height = 9/1.2, width = 16/1.2)       
