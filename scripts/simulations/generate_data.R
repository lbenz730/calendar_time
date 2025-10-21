library(tidyverse)
library(glue)
library(arrow)
library(splines)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Data
df_inform <- read_parquet(glue('{data_dir}/tv_effects/weight_sims_inform.parquet'))


### Function to generate simulated dataset
###
### n_subjects = # of subjects (n)
### n_trials = # of trials (M) 
### covariate_shift = logical, if the population changes over time 
### covariate_shift_fx = list of functions of how the covariates change over time
### outcome_model = outcome model covariate list
### treatment_model = treatment model covariate list 

generate_data <- function(params) {
  ### Unpack simulation parameters 
  n_subjects <- params$n_subjects
  n_trials <- params$n_trials
  covariate_shift <- params$covariate_shift
  covariate_shift_fx <- params$covariate_shift_fx
  treatment_model <- params$treatment_model
  outcome_model <- params$outcome_model
  sigma_rss <- params$sigma_rss
  
  ### (1) Covariates for large sample of eligible subjects we could possible draw from
  subject_ix <- sample(1:nrow(df_inform), n_subjects)
  df_covariates <- 
    df_inform %>% 
    slice(subject_ix) %>% 
    mutate('subject_id' = 1:nrow(.)) %>% 
    select(subject_id, eligible, race, gender, site, baseline_age, baseline_bmi,
           t2dm, insulin, hypertension, hypertension_rx, 
           dyslipidemia, antilipemic_rx, smoking_status)  
  
  ### (2) Expand to Subject Trials
  df_trials <- 
    df_covariates %>% 
    uncount(n_trials) %>% 
    group_by(subject_id) %>% 
    mutate('trial_id' = 1:n_trials) %>% 
    ungroup() %>% 
    select(subject_id, trial_id,  eligible, everything()) %>% 
    inner_join(
      as_tibble(ns(1:n_trials, df = 3)) %>% 
        set_names(paste0('spline_', 1:ncol(.))) %>% 
        mutate('trial_id' = 1:n_trials),
      by = 'trial_id'
    ) 
  
  ### (3) Apply Covariate Shift if Applicable
  if(covariate_shift) {
    for(i in 1:length(covariate_shift_fx)) {
      
      ### Shift Covariate
      pkg <- covariate_shift_fx[[i]]
      var_name <- names(covariate_shift_fx)[i]
      X <- cbind(df_trials[[var_name]], kronecker(matrix(1, n_subjects), pkg$X))
      df_trials[[var_name]] <- as.vector(X %*% pkg$beta)
      
      ### Some covariates the shift is in the mean trend rather than the covariate itself
      if(pkg$random_walk) {
        df_trials[[var_name]] <- df_trials[[var_name]] + rnorm(n = nrow(df_trials), mean = 0, sd = pkg$walk_sd)
      } else if(pkg$probability) { ### Others is the P(covariate = 1) that's trending
        df_trials[[var_name]] <- rbinom(n = nrow(df_trials), size = 1, prob = df_trials[[var_name]])
      }
    }
  }
  
  ### (4) Treatment Assignment
  df_trials <- 
    df_trials %>% 
    mutate('p_treatment' = expit(compute_model(df_trials, treatment_model))) %>% 
    mutate('surgery' = rbinom(n = nrow(df_trials), size = 1, prob = p_treatment)) 
  
  ### (5) Outcome
  df_trials <- 
    df_trials %>% 
    mutate('pct_weight_change' = rnorm(n = nrow(df_trials),
                                       mean = compute_model(df_trials, outcome_model),
                                       sd = sigma_rss))
  
  
  return(df_trials)
}


compute_truth <- function(params) {
  ### Unpack simulation parameters 
  n_subjects <- 2500000
  n_trials <- params$n_trials
  covariate_shift <- params$covariate_shift
  covariate_shift_fx <- params$covariate_shift_fx
  treatment_model <- params$treatment_model
  outcome_model <- params$outcome_model
  sigma_rss <- params$sigma_rss
  
  ### (1) Covariates for large sample of eligible subjects we could possible draw from
  subject_ix <- sample(1:nrow(df_inform), n_subjects)
  df_covariates <- 
    df_inform %>% 
    slice(subject_ix) %>% 
    mutate('subject_id' = 1:nrow(.)) %>% 
    select(subject_id, race, gender, site, baseline_age, baseline_bmi,
           t2dm, insulin, hypertension, hypertension_rx, 
           dyslipidemia, antilipemic_rx, smoking_status)  
  
  ### (2) Expand to Subject Trials
  df_trials <- 
    df_covariates %>% 
    uncount(n_trials) %>% 
    group_by(subject_id) %>% 
    mutate('trial_id' = 1:n_trials) %>% 
    ungroup() %>% 
    select(subject_id, trial_id, everything()) %>% 
    inner_join(
      as_tibble(ns(1:n_trials, df = 3)) %>% 
        set_names(paste0('spline_', 1:ncol(.))) %>% 
        mutate('trial_id' = 1:n_trials),
      by = 'trial_id'
    ) 
  
  
  ### (3) Apply Covariate Shift if Applicable
  if(covariate_shift) {
    for(i in 1:length(covariate_shift_fx)) {
      
      ### Shift Covariate
      pkg <- covariate_shift_fx[[i]]
      var_name <- names(covariate_shift_fx)[i]
      X <- cbind(df_trials[[var_name]], kronecker(matrix(1, n_subjects), pkg$X))
      df_trials[[var_name]] <- as.vector(X %*% pkg$beta)
      
      ### Some covariates the shift is in the mean trend rather than the covariate itself
      if(pkg$random_walk) {
        df_trials[[var_name]] <- df_trials[[var_name]] + rnorm(n = nrow(df_trials), mean = 0, sd = pkg$walk_sd)
      } else if(pkg$probability) { ### Others is the P(covariate = 1) that's trending
        df_trials[[var_name]] <- rbinom(n = nrow(df_trials), size = 1, prob = df_trials[[var_name]])
      }
    }
  }
  
  ### (4) Treatment Assignment
  df_trials <- 
    df_trials %>% 
    mutate('p_treatment' = expit(compute_model(df_trials, treatment_model))) %>% 
    mutate('surgery' = rbinom(n = nrow(df_trials), size = 1, prob = p_treatment)) 
  
  ### (5) Outcome
  df_trials <- 
    df_trials %>% 
    mutate('mu1' = compute_model(mutate(df_trials, 'surgery' = 1), outcome_model),
           'mu0' = compute_model(mutate(df_trials, 'surgery' = 0), outcome_model))
  
  ### Testing it worked
  df_truth <- 
    df_trials %>% 
    group_by(trial_id) %>% 
    summarise('true_ate' = mean(mu1 - mu0)) %>% 
    mutate('trt_change' = params$trt_change,
           'shift' = params$shift)
  
  return(df_truth)
  
}
