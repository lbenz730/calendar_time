library(tidyverse)
library(glue)
library(arrow)
library(splines)
library(mgcv)
library(nleqslv)
library(ranger)
source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R')

set.seed(81923)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
model_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets/tv_effects/weight_models'

### Parameters
args <- commandArgs(trailingOnly = T)
baseline_trial_id <- as.numeric(args[1])

### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 

df_index <- 
  crossing('subject_id' = unique(df_trials$subject_id),
           'trial_id' = unique(df_trials$trial_id))

fit_EIF_chi <- function(df_trials, outcome, n_splits, procedures = c('RYGB', 'SLEEVE', 'AGB'), trial_range = 1:84) {
  
  cat('Computing Empirical Eligibility & Creating Baseline Population for Standardization\n')
  ### Empirical Eligibility
  df_trials$pct_weight_change <- df_trials[[outcome]]
  
  ### Data Set for training models
  df_trials_elig <- 
    df_trials %>% 
    filter(eligible) %>% 
    filter(surgery == 0 | bs_type %in% procedures) %>% 
    filter(!is.na(pct_weight_change)) %>% 
    filter(trial_id %in% trial_range) %>% 
    mutate('site' = as.factor(site),
           'smoking_status' = as.factor(smoking_status),
           'gender' = as.factor(gender),
           'race' = as.factor(race))
  
  ### Baseline Trial Expansion to include entire trial range
  baseline_trials <- 
    df_trials_elig %>% 
    filter(trial_id == baseline_trial_id) %>% 
    group_by(subject_id) %>% 
    reframe('trial_id' = trial_range) %>% 
    inner_join(
      df_trials_elig %>% 
        filter(trial_id == baseline_trial_id) %>% 
        select(-trial_id), 
      by = 'subject_id'
    ) %>% 
    mutate('calendar_year' = 2004 + ceiling(trial_id/12))
  
  ### Empirical Eligibility in Baseline Target Population
  df_elig <- 
    df_index %>% 
    left_join(baseline_trials %>% 
                mutate('eligible' = ifelse(surgery == 1 & !(bs_type %in% procedures), F, eligible)) %>% 
                select(subject_id, trial_id, eligible),
              by = c('subject_id', 'trial_id')) %>% 
    mutate('eligible' = ifelse(is.na(eligible), F, eligible)) %>% 
    group_by(trial_id) %>% 
    summarise('p_elig' = mean(eligible)) %>% 
    filter(trial_id %in% trial_range)
  
  ### Sample splitting
  ### Split by subject to ensure all of a subject's trials end up in the same split\
  cat('Sample Splitting\n')
  subj_ids <- unique(df_trials_elig$subject_id)
  split_ids <- 
    split(sample( 1:length(subj_ids) ), ceiling(1:length(subj_ids) / ceiling(length(subj_ids)/n_splits)))
  
  ### Create Space for Out of Sample Fitted Values
  baseline_trials <- 
    baseline_trials %>% 
    mutate('fold_id' = NA_real_,
           'row_ix' = 1:nrow(.),
           'mu0_hat' = NA_real_,
           'mu1_hat' = NA_real_,
           'pi_hat' = NA_real_,
           'mu0j_hat' = NA_real_,
           'mu1j_hat' = NA_real_,
           'xi_j' = NA_real_,
           'xi_m' = NA_real_)
  
  
  df_trials_elig <- 
    df_trials_elig %>% 
    mutate('fold_id' = NA_real_,
           'row_ix' = 1:nrow(.),
           'pi_hat' = NA_real_)
  
  xi_ratios <- NULL
  for(i in 1:n_splits) {
    test_ids <- subj_ids[ split_ids[[i]] ]
    train_ids <- subj_ids[ setdiff(1:length(subj_ids), split_ids[[i]]) ]
    
    df_train <- 
      df_trials_elig %>% 
      filter(subject_id %in% train_ids)
    
    ### Retain pi_hat on all test for consistent weight truncation
    test_large <- 
      df_trials_elig %>% 
      filter(subject_id %in% test_ids)
    
    ### Hold out predictions for baseline trial population
    df_test <- 
      baseline_trials %>% 
      filter(subject_id %in% test_ids)
    
    ### Load in Models
    cat('Load Models\n')
    proc <- paste(sort(procedures), collapse = '-')
    rf_mu_0 <- read_rds(glue('{model_dir}/{outcome}/{proc}/rf_mu0_holdout_{i}.rds'))
    rf_mu_1 <- read_rds(glue('{model_dir}/{outcome}/{proc}/rf_mu1_holdout_{i}.rds'))
    rf_pi <- read_rds(glue('{model_dir}/{outcome}/{proc}/rf_pi_holdout_{i}.rds'))
    rf_xi <- read_rds(glue('{model_dir}/{outcome}/{proc}/rf_xi_holdout_{i}.rds'))
    
    ### Save Predictions
    cat('Predicting mu0, mu1, pi models on test set [', i, '/', n_splits, ']\n', sep = '')
    df_test$mu0_hat <- predict(rf_mu_0, data = df_test)$predictions
    df_test$mu1_hat <- predict(rf_mu_1, data = df_test)$predictions
    df_test$mu0j_hat <- predict(rf_mu_0, data = df_test %>% mutate('trial_id' = baseline_trial_id))$predictions
    df_test$mu1j_hat <- predict(rf_mu_1, data = df_test %>% mutate('trial_id' = baseline_trial_id))$predictions
    df_test$pi_hat <- predict(rf_pi, data = df_test)$predictions[,2-df_train$surgery[1]]
    
    ### Build out Transport Matrix 
    cat('Build out Transport Matrix test set [', i, '/', n_splits, ']\n', sep = '')
    transport_matrix <- predict(rf_xi, data = test_large)$predictions
    colnames(transport_matrix) <- paste0('xi_', trial_range)
    transport_matrix <- 
      test_large %>% 
      select(subject_id, trial_id) %>% 
      bind_cols(as_tibble(transport_matrix))
    
    transport_df <- 
      transport_matrix %>% 
      group_by(subject_id, trial_id) %>% 
      pivot_longer(cols = contains('xi_'),
                   names_to = 'baseline_trial', 
                   values_to = 'xi_j',
                   names_prefix = 'xi_') %>% 
      mutate('baseline_trial' = as.numeric(baseline_trial)) %>% 
      group_by(subject_id, trial_id) %>% 
      mutate('xi_m' = xi_j[trial_id == baseline_trial]) %>% 
      ungroup()
    
    ### Get everyone's covariate concordance with all treatment times (m) for transporting to this 
    ### population (j)
    transport_df <- 
      transport_df %>% 
      filter(baseline_trial == baseline_trial_id) %>% 
      select(-baseline_trial) 
    
    xi_ratios <- 
      xi_ratios %>% 
      bind_rows(transport_df)
    
    baseline_trials$mu0_hat[df_test$row_ix] <- df_test$mu0_hat 
    baseline_trials$mu1_hat[df_test$row_ix] <- df_test$mu1_hat 
    baseline_trials$mu0j_hat[df_test$row_ix] <- df_test$mu0j_hat 
    baseline_trials$mu1j_hat[df_test$row_ix] <- df_test$mu1j_hat 
    baseline_trials$pi_hat[df_test$row_ix] <- df_test$pi_hat 
    baseline_trials$fold_id[df_test$row_ix] <- i
    
    df_trials_elig$pi_hat[test_large$row_ix] <- predict(rf_pi, data = test_large)$predictions[,2-df_train$surgery[1]]
    df_trials_elig$fold_id[test_large$row_ix] <- i
  }
  
  ### Gather Nuisance functions
  weight_max <- 
    df_trials_elig %>% 
    ### Weights/stabilization
    group_by(fold_id) %>% 
    mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-pi_hat),
                                     surgery == 1 ~ 1/pi_hat),
           'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat),
                                      surgery == 1 ~ mean(surgery)/(pi_hat))) %>% 
    ungroup() %>% 
    ### Trim Weights at 99th quantile
    group_by(surgery, fold_id) %>%
    summarise('w_treatment_bound' = max(winsorize(w_treatment, q = c(0, 0.99))),
              'sw_treatment_bound' = max(winsorize(sw_treatment, q = c(0, 0.99)))) %>%
    ungroup()
  
  baseline_trials <- 
    baseline_trials %>% 
    mutate('muA_hat' = ifelse(surgery == 1, mu1_hat, mu0_hat),
           'muAj_hat' = ifelse(surgery == 1, mu1j_hat, mu0j_hat)) %>% 
    ### Weights/stabilization
    group_by(surgery, fold_id) %>%
    mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-pi_hat),
                                     surgery == 1 ~ 1/pi_hat),
           'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat),
                                      surgery == 1 ~ mean(surgery)/(pi_hat))) %>% 
    inner_join(weight_max, by = c('surgery', 'fold_id')) %>% 
    ### Trim Weights at 99th quantile (overall)
    mutate('w_treatment' = pmin(w_treatment, w_treatment_bound),
           'sw_treatment' = pmin(sw_treatment, sw_treatment_bound)) %>% 
    ungroup()
  
  ### Read in saved out fold ids for this analysis to assign never eligible to the same fold
  fold_ids <- 
    read_parquet(glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome}_{paste(sort(procedures), collapse = "-")}.parquet')) %>% 
    select(subject_id, trial_id, fold_id)
  
  ### Trial Specific Chi Hat
  baseline_trials <- 
    baseline_trials %>% 
    mutate('chi_IF' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muA_hat),
           'chi_IFj' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muAj_hat))
  
  df_chi <- 
    baseline_trials %>% 
    group_by(trial_id) %>%
    summarise('chi_DR' = mean(chi_IF),
              'chi_DRj' = mean(chi_IFj),
              'chi_gformula' = mean(mu1_hat - mu0_hat),
              'p_treatment' = mean(surgery))
  
  ### Note chi_IF and chi_IFj aren't eligibility weighted
  df_IF <- 
    baseline_trials %>% 
    select(subject_id, trial_id, eligible, 
           surgery, mu1_hat, 
           mu0_hat, pi_hat, chi_IF, mu1j_hat, mu0j_hat, chi_IFj) %>% 
    right_join(df_index %>% filter(trial_id %in% trial_range), by = c('subject_id', 'trial_id')) %>% 
    left_join(fold_ids, by = c('subject_id', 'trial_id')) %>% 
    mutate('eligible' = as.numeric(replace(eligible, is.na(eligible), F)),
           'chi_IF' = replace(chi_IF, is.na(chi_IF), 0),
           'chi_IFj' = replace(chi_IFj, is.na(chi_IFj), 0))
  
  ### IF for cross parameters
  ### Gather eligibility, mu, pi and outcome/treatment for all m
  IF_m <- read_parquet(glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome}_{paste(sort(procedures), collapse = "-")}.parquet'))
  IF_m <- 
    IF_m %>% 
    select(subject_id, trial_id, 'eligible_m' = eligible, 'surgery_m' = surgery,
           'mu1m_hat' = mu1_hat, 'mu0m_hat' = mu0_hat, 'pim_hat' = pi_hat) %>% 
    left_join(
      df_trials %>% 
        select(subject_id, trial_id, 'pct_weight_change_m' = pct_weight_change),
      by = c('subject_id', 'trial_id') 
    )
  
  df_IF_cross <- 
    baseline_trials %>% 
    select(subject_id, trial_id, 'eligible_j' = eligible, surgery, mu1_hat, mu0_hat) %>% 
    right_join(df_index %>% filter(trial_id %in% trial_range), by = c('subject_id', 'trial_id')) %>% 
    left_join(fold_ids, by = c('subject_id', 'trial_id')) %>% 
    inner_join(IF_m, by = c('subject_id', 'trial_id')) %>% 
    mutate('eligible_j' = as.numeric(replace(eligible_j, is.na(eligible_j), F))) %>% 
    group_by(trial_id) %>% 
    mutate('p_elig_j' = mean(eligible_j),
           'p_elig_m' = mean(eligible_m)) %>% 
    ungroup() %>% 
    mutate('muAm_hat' = ifelse(surgery_m == 1, mu1m_hat, mu0m_hat),
           'w_treatment_m' = case_when(surgery_m == 0 ~ 1/(1-pim_hat),
                                       surgery_m == 1 ~ 1/pim_hat),
           'sw_treatment_m' = case_when(surgery_m == 0 ~ mean(1-surgery_m)/(1-pim_hat),
                                        surgery_m == 1 ~ mean(surgery_m)/(pim_hat))) %>% 
    left_join(weight_max, by = c('surgery_m' = 'surgery', 'fold_id')) %>% 
    ### Trim Weights at 99th quantile (overall)
    mutate('w_treatment_m' = pmin(w_treatment_m, w_treatment_bound),
           'sw_treatment_m' = pmin(sw_treatment_m, sw_treatment_bound)) %>% 
    left_join(xi_ratios, by = c('subject_id', 'trial_id')) %>% 
    mutate('xi_ratio' = xi_j/xi_m) %>% 
    group_by(trial_id) %>% 
    mutate('xi_ratio' = winsorize(xi_ratio, q = c(0, 0.99))) %>% 
    ungroup() %>% 
    mutate('IF_j' = eligible_j/p_elig_j * (mu1_hat - mu0_hat),
           'IF_j' = replace(IF_j, eligible_j == 0, 0),
           'IF_m' = eligible_m/p_elig_m * (2 * surgery_m - 1) * w_treatment_m * xi_ratio * p_elig_m/p_elig_j * (pct_weight_change_m - muAm_hat),
           'IF_m' = replace(IF_m, eligible_m == 0, 0),
           'IF_cross' = IF_j + IF_m)
  
  chi_cross <- 
    df_IF_cross %>% 
    group_by(trial_id) %>% 
    summarise('chi_cross' = mean(IF_cross))
  
  df_info <- 
    df_chi %>% 
    inner_join(chi_cross, by = 'trial_id') %>% 
    inner_join(df_elig, by = 'trial_id') %>% 
    select(trial_id, contains('chi'), contains('p_')) %>% 
    mutate('outcome' = outcome,
           'procedures' = paste(sort(procedures), collapse = '/'))
  
  
  
  return( list('df_info' = df_info,
               'df_IF_cross' = df_IF_cross,
               'df_IF' = df_IF) )
}


### Run analysis for each outcome/procedure(s) combination
df_results <- NULL
df_importance <- NULL
df_seeds <- 
  crossing('outcome' = c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr'),
           'procedures' = c('ALL', 'RYGB', 'SLEEVE')) %>% 
  mutate('seed' = sample(1:100000000, 12))

procedure_list <- c('ALL', 'RYGB', 'SLEEVE')
if(baseline_trial_id <= 29) {
  procedure_list <- c('ALL', 'RYGB')
}

for(outcome in c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr')) {
  cat('\nAnalysis for Outcome:', outcome, '\n')
  for(procedures in procedure_list) {
    set.seed(df_seeds$seed[df_seeds$outcome == outcome & df_seeds$procedures == procedures])
    cat('\nAnalysis for Procedure(s):', procedures, '\n')  
    t_range <- 1:84
    
    if(procedures == 'ALL') {
      procedures <-  c('AGB', 'RYGB', 'SLEEVE')
    } else if(procedures == 'SLEEVE') {
      t_range <- 30:84 
    }
    
    ### Run Analysis
    pkg <- 
      fit_EIF_chi(df_trials = df_trials, 
                  outcome = outcome,
                  n_splits = 2, 
                  procedures = procedures,
                  trial_range = t_range)
    
    ### Save Results
    df_results <- bind_rows(df_results, pkg$df_info)
    
    if(!dir.exists(glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}'))) {
      dir.create(glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}'))
    }
    
    write_parquet(pkg$df_IF, glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}/subject_IF_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
    write_parquet(pkg$df_IF_cross, glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}/subject_IF_cross_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
    
    cat('\n\n---------------------------------\n\n')
  }
}

### Can't do Sleeve/RYGB comparison on a population before Sleeve really existed in our data
if(baseline_trial_id >= 30) {
  ### RYGB vs. Sleeve
  ### For ease of analysis pipeline here we let surgery = 1 denote RYGB and 
  ### surgery = 0 denote SLEEVE
  set.seed(10293)
  df_surg <- 
    df_trials %>% 
    mutate('surgery_original' = surgery) %>% 
    mutate('eligible' = ifelse(surgery_original == 0, F, eligible)) %>% 
    mutate('surgery' = case_when(bs_type == 'CONTROL' ~ NA,
                                 bs_type == 'AGB' ~ NA,
                                 bs_type == 'RYGB' ~ 1,
                                 bs_type == 'SLEEVE' ~ 0))
  
  
  
  for(outcome in c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr')) {
    cat('\nAnalysis for Outcome:', outcome, '[RYGB vs. SLEEVE] \n')
    
    
    ### Run Analysis
    pkg <- 
      fit_EIF_chi(df_trials = df_surg, 
                  outcome = outcome,
                  n_splits = 2, 
                  procedures = c('RYGB', 'SLEEVE'),
                  trial_range = 30:84)
    
    ### Save Results
    df_results <- bind_rows(df_results, pkg$df_info)
    
    if(!dir.exists(glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}'))) {
      dir.create(glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}'))
    }
    
    write_parquet(pkg$df_IF, glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}/subject_IF_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
    write_parquet(pkg$df_IF_cross, glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}/subject_IF_cross_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
    
    cat('\n\n---------------------------------\n\n')
  }
}

### Save Results
write_csv(df_results, glue('{data_dir}/tv_effects/weight_DR_results_baseline_trial_{baseline_trial_id}.csv'))
