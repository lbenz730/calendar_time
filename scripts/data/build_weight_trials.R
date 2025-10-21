library(tidyverse)
library(lubridate)
library(arrow)
library(glue)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Demographics, labs etc.
df_subjects <- read_parquet(glue('{data_dir}/microvascular_tte/subjects.parquet'))
weights <- read_parquet(glue('{data_dir}/all_weights_further_cleaned.parquet')) ### Cleaned Weights w/ outliers removed
df_pregnancy <- read_parquet(glue('{data_dir}/microvascular_tte/pregnancy.parquet'))
df_enrollment <- read_parquet(glue('{data_dir}/microvascular_tte/enrollment.parquet'))
diabetes_rx <- read_parquet(glue('{ehr_dir}/parquet_files/combined/diabetes_rx.parquet'))
diabetes_labs <- read_parquet(glue('{data_dir}/microvascular_tte/diabetes_labs.parquet'))
diabetes_dx <- read_parquet(glue('{ehr_dir}/parquet_files/combined/diabetes_dx.parquet'))
smoking <- read_parquet(glue('{data_dir}/microvascular_tte/smoking.parquet'))
surgical_px <- read_parquet( glue('{data_dir}/bs_types_reviewed.parquet'))
obesity_comorbidities <- read_parquet(glue('{data_dir}/bariatric_tte/obesity_comorbidities.parquet'))
antilipemic_rx <- read_parquet(glue('{data_dir}/antilipemics_rx.parquet'))
hypertension_rx <- read_parquet(glue('{data_dir}/hypertension_rx.parquet'))

### Drop measures on the same date
weights <-
  weights %>% 
  distinct(subject_id, measure_date, .keep_all = T) %>% 
  filter(!is.na(height)) %>% 
  filter(bmi >= 10) ## Extreme Outliers

smoking <- 
  smoking %>% 
  distinct(subject_id, contact_date, .keep_all = T)

diabetes_labs <- 
  diabetes_labs %>% 
  distinct(subject_id, lab_date, .keep_all = T)

### Function to build a trial (pre-discrete time)
### trial_id = which trial in sequence to build
### study_start = start date of trial sequence
### study_end = end date of study (administrative end of study)
### lookbacks = list of lookback window for various aspects of the studies
###     a1c = # of months to lookback to ascertain T2DM
###     bmi = # of months to lookback for BMI
###     t2dm_rx = # of months to lookback for T2DM rx
###     t2dm_icd9 = # of months to lookback for T2DM via ICD-9 among patients on METFORMIN as the sole indicator of diabetes
###     htdys = # of months to lookback for hypertension + dyslipidemia dx
###     egfr = # of months to look back for estimated glucose filtration rate
build_trial <- function(trial_id, study_start, study_end, lookbacks) {
  ### Trial 'Enrollment' Period
  trial_start <- study_start %m+% months(trial_id - 1)
  trial_end <- study_start %m+% months(trial_id) - 1
  
  ### Pregnancies to exclude in the last year
  pregancy_1yr <- 
    df_pregnancy %>% 
    filter(adate <= trial_start,
           adate >= trial_start %m-% years(1)) %>% 
    distinct(subject_id)
  
  ### BS Types of Patients Getting Surgery this month
  bs_types <- 
    surgical_px %>% 
    filter(index_date >= trial_start & index_date <= trial_end) %>% 
    select(subject_id, bs_type)
  
  
  trial_population <-
    df_subjects %>% 
    
    ### Remove pregnancies in last year
    anti_join(pregancy_1yr, by = 'subject_id') %>% 
    
    ### Remove old cases and code surgery in month
    filter(is.na(index_date) | index_date >= trial_start) %>% 
    mutate('surgery' = as.numeric(!is.na(index_date) & index_date >= trial_start & index_date <= trial_end)) %>% 
    
    ### Remove surgeries that aren't primary BS
    left_join(bs_types, by = 'subject_id') %>% 
    mutate('bs_type' = ifelse(surgery == 0, 'CONTROL', bs_type)) %>% 
    filter(bs_type != 'EXCLUDE') %>% 
    
    ### Enrollment (Must be enrolled for at least 1 year before trial start)
    inner_join(df_enrollment %>% 
                 filter(enr_1yr <= trial_start, enr_end >= trial_start),
               by = 'subject_id') 
  
  ### Determine Eligibility
  ### BMI
  ### max(BMI) >= 35 in last year (surgery) and 34+ at date of surgery
  ### most recent BMI >= 35 at most recent year (control) 
  ###
  df_bmi <- 
    weights %>% 
    filter(measure_date >= trial_start %m-% months(lookbacks$bmi),
           measure_date <= pmin(index_date, trial_end, na.rm = T)) %>% 
    arrange(measure_date) %>% 
    group_by(subject_id) %>% 
    summarise('baseline_bmi' = last(bmi),
              'max_bmi' = max(bmi))
  
  df_a1c <- 
    diabetes_labs %>% 
    left_join(surgical_px, by = 'subject_id') %>%
    filter(lab_date >= trial_start %m-% months(lookbacks$a1c),
           lab_date <= pmin(index_date, trial_end, na.rm = T)) %>% 
    ### Remove extreme A1c Values
    filter((test_type == 'HGBA1C' & result <= 14) | (test_type == 'GLU_F' & (result + 46.7)/28.7 <= 14)) %>% 
    arrange(lab_date) %>% 
    group_by(subject_id) %>% 
    summarise('baseline_a1c' = last(result[test_type == 'HGBA1C']),
              'baseline_gluF' = last(result[test_type == 'GLU_F'])) %>% 
    mutate('baseline_a1c' = ifelse(is.na(baseline_a1c), (baseline_gluF + 46.7)/28.7, baseline_a1c),
           'baseline_gluF' = ifelse(is.na(baseline_gluF), baseline_a1c * 28.7 - 46.7, baseline_gluF))
  
  ### Active prescriptions
  df_rx <- 
    diabetes_rx %>% 
    select(-index_date) %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(rxdate <= pmin(index_date, trial_end, na.rm = T),
           rxdate + rxsup >= trial_start %m-% months(lookbacks$rx)) %>% 
    group_by(subject_id) %>%
    summarise('t2dm_rx' = 1,
              'insulin' = max(0, insulin_flg, na.rm = T),
              'metformin_only' = as.numeric(all(generic == 'METFORMIN')))
  
  ### Diabetes Diagnosis for subjects w/ metformin as the sole indicator of T2DM
  df_dx <-
    diabetes_dx %>% 
    select(-index_date) %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(adate <= pmin(index_date, trial_end, na.rm = T),
           adate >= trial_start %m-% months(lookbacks$t2dm_icd9)) %>% 
    filter(grepl('250', dx)) %>% 
    group_by(subject_id) %>% 
    summarise('diabetes_icd9' = 1) %>% 
    ungroup()
  
  ### Most recent Smoking Status
  df_smoking <- 
    smoking %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(contact_date <= pmin(trial_end, index_date, na.rm = T)) %>% 
    arrange(desc(contact_date)) %>% 
    group_by(subject_id)  %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(-index_date, -bs_type)
  
  ### Processing of Covariates
  df_obesity <- 
    obesity_comorbidities %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(adate <= pmin(index_date, trial_end, na.rm = T),
           adate >= trial_start %m-% months(lookbacks$htdys)) %>% 
    group_by(subject_id) %>% 
    summarise('hypertension' = max(hypertension),
              'dyslipidemia' = max(dyslipidemia))
  
  
  ### Hypertension + Antilipemic Rx
  ht_rx <- 
    hypertension_rx %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(rxdate <= pmin(index_date, trial_end, na.rm = T),
           rxdate + rxsup >= trial_start %m-% months(lookbacks$rx)) %>% 
    group_by(subject_id) %>%
    summarise('hypertension_rx' = 1)
  
  dys_rx <- 
    antilipemic_rx %>% 
    left_join(surgical_px, by = 'subject_id') %>% 
    filter(rxdate <= pmin(index_date, trial_end, na.rm = T),
           rxdate + rxsup >= trial_start %m-% months(lookbacks$rx)) %>% 
    group_by(subject_id) %>%
    summarise('antilipemic_rx' = 1)
  
  ### Dataset for trial (before expansion)
  df_trial <- 
    trial_population %>% 
    mutate('trial_id' = trial_id,
           'trial_start' = trial_start,
           'calendar_year' = year(trial_start)) %>% 
    mutate('baseline_age' = as.numeric(trial_start - birth_date)/365.25) %>% 
    left_join(df_bmi, by = 'subject_id') %>% 
    left_join(df_a1c, by = 'subject_id') %>% 
    left_join(df_rx, by = 'subject_id') %>% 
    left_join(df_dx, by = 'subject_id') %>% 
    left_join(df_smoking, by = 'subject_id') %>% 
    left_join(df_obesity, by = 'subject_id') %>% 
    left_join(ht_rx, by = 'subject_id') %>% 
    left_join(dys_rx, by = 'subject_id') %>% 
    ### Determine T2DM Status
    mutate('t2dm' = case_when(baseline_a1c >= 6.5 | baseline_gluF >= 126 ~ 1, ### Via Labs
                              t2dm_rx == 1 & !metformin_only ~ 1, ### Via Rx other than metformin
                              t2dm_rx == 1 & metformin_only & diabetes_icd9 == 1 ~ 1, ### Via Rx (metformin) + 250.x ICD-9
                              is.na(t2dm_rx) & baseline_a1c < 6.5 & baseline_gluF < 126 ~ 0,
                              T ~ NA)) %>% 
    ### Eligibility Status 
    mutate('eligible' = (surgery == 1 & max_bmi >= 35 & baseline_bmi >= 33) | (surgery == 0 & baseline_bmi >= 35) ) 
  
  ### Light Data Cleaning
  df_trial <-
    df_trial %>% 
    mutate('race' = ifelse(race %in% c('MU', 'OT', 'UN', 'AS', 'IN', 'HP'), 'OT', race)) %>% 
    mutate('dyslipidemia' = ifelse(is.na(dyslipidemia), 0, dyslipidemia),
           'hypertension' = ifelse(is.na(hypertension), 0, hypertension),
           't2dm_rx' = ifelse(is.na(t2dm_rx), 0, t2dm_rx),
           'antilipemic_rx' = ifelse(is.na(antilipemic_rx), 0, antilipemic_rx),
           'hypertension_rx' = ifelse(is.na(hypertension_rx), 0, hypertension_rx),
           'insulin' = ifelse(is.na(insulin), 0, insulin),
           'smoking_status' = ifelse(is.na(smoking_status), 'no_self_report', smoking_status)) %>% 
    ### Remove Subjects w/out Baseline BMI and only keep those w/ known T2DM Status 
    filter(!is.na(t2dm),
           !is.na(baseline_bmi))
  
  ### 6 month Weight (+/- 2 months)
  bmi_6mo <- 
    weights %>% 
    filter(measure_date >= trial_start %m+% months(4),
           measure_date <= trial_start %m+% months(8)) %>% 
    filter(subject_id %in% df_trial$subject_id[df_trial$eligible]) %>% 
    mutate('mo6_date' = trial_start %m+% months(6)) %>% 
    mutate('time' = as.numeric(measure_date - trial_start),
           'time_6mo' = as.numeric(mo6_date - trial_start)) %>% 
    group_by(subject_id, index_date, time_6mo) %>% 
    nest() %>% 
    mutate('bmi_6mo' = map2_dbl(data, time_6mo, ~bmi_assess(.x, .y))) %>% 
    ungroup() %>%
    select(subject_id, bmi_6mo)  
  
  ### 1 Year Weight (+/- 3 months)
  bmi_1yr <- 
    weights %>% 
    filter(measure_date >= trial_start %m+% months(9),
           measure_date <= trial_start %m+% months(15)) %>% 
    filter(subject_id %in% df_trial$subject_id[df_trial$eligible]) %>% ### Only compute outcomes for eligible subjects
    mutate('yr1_date' = trial_start %m+% months(12)) %>% 
    mutate('time' = as.numeric(measure_date - trial_start),
           'time_1yr' = as.numeric(yr1_date - trial_start)) %>% 
    group_by(subject_id, index_date, time_1yr) %>% 
    nest() %>% 
    mutate('bmi_1yr' = map2_dbl(data, time_1yr, ~bmi_assess(.x, .y))) %>% 
    ungroup() %>%
    select(subject_id, bmi_1yr)  
  
  
  ### 2 Year Weight (+/- 6 months)
  bmi_2yr <- 
    weights %>% 
    filter(measure_date >= trial_start %m+% months(18),
           measure_date <= trial_start %m+% months(30)) %>% 
    filter(subject_id %in% df_trial$subject_id[df_trial$eligible]) %>% 
    mutate('yr2_date' = trial_start %m+% months(24)) %>% 
    mutate('time' = as.numeric(measure_date - trial_start),
           'time_2yr' = as.numeric(yr2_date - trial_start)) %>% 
    group_by(subject_id, index_date, time_2yr) %>% 
    nest() %>% 
    mutate('bmi_2yr' = map2_dbl(data, time_2yr, ~bmi_assess(.x, .y))) %>% 
    ungroup() %>%
    select(subject_id, bmi_2yr)  
  
  ### 3 Year Weight (+/- 6 months)
  bmi_3yr <- 
    weights %>% 
    filter(measure_date >= trial_start %m+% months(30),
           measure_date <= trial_start %m+% months(42)) %>% 
    filter(subject_id %in% df_trial$subject_id[df_trial$eligible]) %>% 
    mutate('yr3_date' = trial_start %m+% months(36)) %>% 
    mutate('time' = as.numeric(measure_date - trial_start),
           'time_3yr' = as.numeric(yr3_date - trial_start)) %>% 
    group_by(subject_id, index_date, time_3yr) %>% 
    nest() %>% 
    mutate('bmi_3yr' = map2_dbl(data, time_3yr, ~bmi_assess(.x, .y))) %>% 
    ungroup() %>%
    select(subject_id, bmi_3yr)
  
  ### Bind all outcome
  df_outcomes <- 
    bmi_6mo %>% 
    full_join(bmi_1yr, by = 'subject_id') %>% 
    full_join(bmi_2yr, by = 'subject_id') %>% 
    full_join(bmi_3yr, by = 'subject_id') 
  
  
  df_trial <- 
    df_trial %>% 
    left_join(df_outcomes, by = c('subject_id')) %>% 
    filter(gender != 'O') %>% 
    mutate('delta_6mo' = (bmi_6mo - baseline_bmi)/baseline_bmi,
           'delta_1yr' = (bmi_1yr - baseline_bmi)/baseline_bmi,
           'delta_2yr' = (bmi_2yr - baseline_bmi)/baseline_bmi,
           'delta_3yr' = (bmi_3yr - baseline_bmi)/baseline_bmi) %>% 
    select(subject_id, trial_id, calendar_year, eligible, ### IDs
           surgery, bs_type, index_date, ### Treatment
           baseline_bmi, bmi_6mo, bmi_1yr, bmi_2yr, bmi_3yr, ### BMI Related Outcomes
           delta_6mo, delta_1yr, delta_2yr, delta_3yr, 
           gender, race, site, baseline_age, 
           t2dm, insulin,  ### T2DM Covariates
           hypertension, hypertension_rx, dyslipidemia, antilipemic_rx, smoking_status ### Other Covariates
    )
  
  return(df_trial)
}

### Process Outcomes (Weight Change)
### Method of (Thaweethai et al., 2021)
bmi_assess <- function(df, time_) {
  ### If single weight in period use that 
  if(nrow(df) == 1) {
    return(df$bmi) 
  } 
  
  ### If all measures occur on one side of X year date use closest
  if(all(df$time <= time_) | all(df$time >= time_)) {
    return( df$bmi[which.min(abs(df$time - time_))] )
  }
  
  ### Otherwise Use Linear Regression to Estimate
  bmi_model <- lm(bmi ~ time, data = df)
  return(predict(bmi_model, newdata = tibble('time' = time_)))
  
}



### Build Trial
args <- commandArgs(trailingOnly = T)
t_id <- as.numeric(args[1])

df_trial <- 
  build_trial(study_start = as.Date('2005-01-01'),
              study_end = as.Date('2015-09-30'),
              lookbacks =   list('a1c' = 24, ### 2 years for A1c
                                 'rx' = 0, ### Active Rx
                                 't2dm_icd9' = 12, 
                                 'bmi' = 12, ### 1 year for BMI
                                 'htdys' = 12),
              trial_id = t_id)


### Save file
if(!dir.exists(glue('{data_dir}/tv_effects/weight_trials'))) {
  dir.create(glue('{data_dir}/tv_effects/weight_trials'))
}
write_parquet(df_trial, glue('{data_dir}/tv_effects/weight_trials/trial_{t_id}.parquet'))
