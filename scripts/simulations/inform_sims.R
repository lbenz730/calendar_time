library(tidyverse)
library(glue)
library(arrow)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 

df_inform <- 
  df_trials %>% 
  filter(eligible == 1) %>% 
  filter(!is.na(delta_3yr)) %>% 
  mutate('pct_weight_change' = delta_3yr) %>% 
  mutate('gender' = ifelse(gender == 'M', 1, 0),
         'race' = ifelse(race == 'WH', 1, 0))

### Fit some Models to get rough coefficient estimates
propensity_model <- 
  speedglm::speedglm(surgery ~ 
                       baseline_bmi + gender + race + site + baseline_age + 
                       t2dm + insulin + hypertension + hypertension_rx + 
                       dyslipidemia + antilipemic_rx + smoking_status,
                     family = binomial(),
                     data = df_inform)

outcome_model <- 
  speedglm::speedlm(pct_weight_change ~
                      surgery + baseline_bmi + gender + race + site + baseline_age + 
                      t2dm + insulin + hypertension + hypertension_rx + 
                      dyslipidemia + antilipemic_rx + smoking_status,
                    data = df_inform)

cat(paste("'", names(coefficients(propensity_model)), "' = ", 
      sprintf('%0.5f', coefficients(propensity_model)), sep = ''), sep = ',\n')

cat(paste("'", names(coefficients(outcome_model)), "' = ", 
          sprintf('%0.5f', coefficients(outcome_model)), sep = ''), sep = ',\n')


df_inform <- 
  df_inform %>% 
  filter(smoking_status != 'no_self_report') 
write_parquet(df_inform, glue('{data_dir}/tv_effects/weight_sims_inform.parquet')) 


