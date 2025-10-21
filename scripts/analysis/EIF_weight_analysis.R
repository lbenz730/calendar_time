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
save_importance <- F

### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 

df_index <- 
  crossing('subject_id' = unique(df_trials$subject_id),
           'trial_id' = unique(df_trials$trial_id))

### Function to estimate the EIF of Chi
###
### df_trials = dataset of all combined target trials
### ipw_formula = formula for pi
### outcome_formula = formula for mu
### transport_formula = formula for xi
### outcome = outcome (among delta_6mo, delta_1yr, delta_2yr, delta_3yr)
### n_splits = # of outer splits for sample splitting 
### n_subsplits = # of inner splits (needed for computing the projection + loss function evaluation)
### procedures = vector of procedures to include in surgical comparison
### importance = T/F, whether or not to save variable importance
### trial_range = which trial(s) to include in analysis
###
### Note the current set-up is designed for 2 outer folds (e.g. so we can 
### save everything and dump the file and do projection analysis after)
### We'd need to alter bookkeeping if more than 2 outer splits were used
### Any number of internal splits can be used for sample splitting in estimation
fit_EIF_chi <- function(df_trials, ipw_formula, outcome_formula, transport_formula, outcome, 
                        n_splits, n_subsplits,
                        procedures = c('RYGB', 'SLEEVE', 'AGB'),
                        importance = T,
                        trial_range = 1:84) {
  
  cat('Computing Empirical Eligibility\n')
  ### Empirical Eligibility
  df_trials$pct_weight_change <- df_trials[[outcome]]
  
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
  
  df_elig <- 
    df_index %>% 
    left_join(df_trials_elig %>% 
                mutate('eligible' = ifelse(surgery == 1 & !(bs_type %in% procedures), F, eligible)) %>% 
                select(subject_id, trial_id, eligible),
              by = c('subject_id', 'trial_id')) %>% 
    mutate('eligible' = ifelse(is.na(eligible), F, eligible)) %>% 
    group_by(trial_id) %>% 
    summarise('p_elig' = mean(eligible)) %>% 
    filter(trial_id %in% trial_range)
  
  ### Sample splitting
  ### Split by subject to ensure all of a subject's trials end up in the same split
  cat('Sample Splitting\n')
  subj_ids <- unique(df_trials_elig$subject_id)
  split_ids <- 
    split(sample( 1:length(subj_ids) ), ceiling(1:length(subj_ids) / ceiling(length(subj_ids)/n_splits)))
  
  ### Create Space for Out of Sample Fitted Values
  df_trials_elig <- 
    df_trials_elig %>% 
    mutate('fold_id' = NA_real_,
           'sub_fold_id' = NA_real_,
           'row_ix' = 1:nrow(.),
           'mu0_hat' = NA_real_,
           'mu1_hat' = NA_real_,
           'pi_hat' = NA_real_,
           'mu0_hat_sub' = NA_real_,
           'mu1_hat_sub' = NA_real_,
           'pi_hat_sub' = NA_real_)
  
  df_importance <- NULL
  
  for(i in 1:n_splits) {
    test_ids <- subj_ids[ split_ids[[i]] ]
    train_ids <- subj_ids[ setdiff(1:length(subj_ids), split_ids[[i]]) ]
    
    df_train <- 
      df_trials_elig %>% 
      filter(subject_id %in% train_ids)
    
    df_test <- 
      df_trials_elig %>% 
      filter(subject_id %in% test_ids)
    
    ### Outer Level Fitting 
    cat('Fitting Outcome Model on split [', i, '/', n_splits, ']\n', sep = '')
    ### Fit Outcome Model
    ### Stratification to deal w/ imbalanced
    
    ### Note we do difference RF for interpretation and estimation because impurity corrected (slightly) modifies the splitting criterion
    rf_mu_0 <- 
      ranger(outcome_formula, 
             num.threads = 64,
             max.depth = 10,
             num.trees = 500,
             importance = 'none',
             data = df_train %>% filter(surgery == 0),
             seed = 109)
    
    rf_mu_1 <- 
      ranger(outcome_formula, 
             num.threads = 64,
             max.depth = 10,
             num.trees = 500,
             importance = 'none',
             data = df_train %>% filter(surgery == 1),
             seed = 11031)
    
    
    ### Fit Propensity Model
    cat('Fitting PS Model on split [', i, '/', n_splits, ']\n', sep = '')
    rf_pi <- 
      ranger(ipw_formula, 
             num.threads = 64,
             max.depth = 2, ### Cap Depth 
             num.trees = 500,
             probability = T,
             importance = 'none',
             data = df_train,
             seed = 81002)
    
    cat('Fitting Transport Model on split [', i, '/', n_splits, ']\n', sep = '')
    rf_xi <- 
      ranger(transport_formula, 
             num.threads = 64,
             max.depth = 10, ### Cap Depth 
             num.trees = 500,
             probability = T,
             importance = 'none',
             data = df_train,
             seed = 819992)
    
    ### Save Predictions and Models
    cat('Predicting mu0, mu1, pi models on test set [', i, '/', n_splits, ']\n', sep = '')
    df_test$mu0_hat <- predict(rf_mu_0, data = df_test)$predictions
    df_test$mu1_hat <- predict(rf_mu_1, data = df_test)$predictions
    df_test$pi_hat <- predict(rf_pi, data = df_test)$predictions[,2-df_train$surgery[1]]
    
    df_trials_elig$mu0_hat[df_test$row_ix] <- df_test$mu0_hat 
    df_trials_elig$mu1_hat[df_test$row_ix] <- df_test$mu1_hat 
    df_trials_elig$pi_hat[df_test$row_ix] <- df_test$pi_hat 
    df_trials_elig$fold_id[df_test$row_ix] <- i ### Save the index of the holdout fold(s)
    
    proc <- paste(sort(procedures), collapse = '-')
    if(!dir.exists(glue('{model_dir}/{outcome}/{proc}'))) {
      dir.create(glue('{model_dir}/{outcome}/{proc}'), recursive = T)
    }
    write_rds(rf_mu_0, glue('{model_dir}/{outcome}/{proc}/rf_mu0_holdout_{i}.rds'))
    write_rds(rf_mu_1, glue('{model_dir}/{outcome}/{proc}/rf_mu1_holdout_{i}.rds'))
    write_rds(rf_pi, glue('{model_dir}/{outcome}/{proc}/rf_pi_holdout_{i}.rds'))
    write_rds(rf_xi, glue('{model_dir}/{outcome}/{proc}/rf_xi_holdout_{i}.rds'))
    
    ### Sub-Split the Training Set to Get Hold out predictions (necessary step for evaluating projections later)
    train_subj_ids <- unique(df_train$subject_id)
    subsplit_ids <- 
      split(sample( 1:length(train_subj_ids) ), ceiling(1:length(train_subj_ids) / ceiling(length(train_subj_ids)/n_subsplits)))
    
    for(j in 1:n_subsplits) {
      sub_test_ids <- train_subj_ids[ subsplit_ids[[j]] ]
      sub_train_ids <- train_subj_ids[ setdiff(1:length(train_subj_ids), subsplit_ids[[j]]) ]
      
      df_train_sub <- 
        df_train %>% 
        filter(subject_id %in% sub_train_ids)
      
      df_test_sub <- 
        df_train %>% 
        filter(subject_id %in% sub_test_ids)
      
      cat('Fitting Outcome Model on subsplit [', j, '/', n_subsplits, '] for training fold ', i, '\n', sep = '')
      ### Note we do difference RF for interpretation and estimation because impurity corrected (slightly) modifies the splitting criterion
      rf_mu_0_sub <- 
        ranger(outcome_formula, 
               num.threads = 64,
               max.depth = 10,
               num.trees = 500,
               importance = 'none',
               data = df_train_sub %>% filter(surgery == 0),
               seed = 12109)
      
      
      rf_mu_1_sub <- 
        ranger(outcome_formula, 
               num.threads = 64,
               max.depth = 10,
               num.trees = 500,
               importance = 'none',
               data = df_train_sub %>% filter(surgery == 1),
               seed = 711031)
      
      
      ### Fit Propensity Model
      cat('Fitting PS Model on subsplit [', j, '/', n_subsplits, '] for training fold ', i, '\n', sep = '')
      rf_pi_sub <- 
        ranger(ipw_formula, 
               num.threads = 64,
               max.depth = 2, ### Cap Depth 
               num.trees = 500,
               probability = T,
               importance = 'none',
               data = df_train_sub,
               seed = 781002)
      
      ### Save Predictions
      cat('Predicting mu0, mu1, pi models on test subsplit [', j, '/', n_subsplits, '] for training fold ', i, '\n', sep = '')
      df_test_sub$mu0_hat <- predict(rf_mu_0_sub, data = df_test_sub)$predictions
      df_test_sub$mu1_hat <- predict(rf_mu_1_sub, data = df_test_sub)$predictions
      df_test_sub$pi_hat <- predict(rf_pi_sub, data = df_test_sub)$predictions[,2-df_train_sub$surgery[1]]
      
      df_trials_elig$mu0_hat_sub[df_test_sub$row_ix] <- df_test_sub$mu0_hat 
      df_trials_elig$mu1_hat_sub[df_test_sub$row_ix] <- df_test_sub$mu1_hat 
      df_trials_elig$pi_hat_sub[df_test_sub$row_ix] <- df_test_sub$pi_hat 
      df_trials_elig$sub_fold_id[df_test_sub$row_ix] <- j ### Save the index of the holdout fold(s)
    }
    
    
    ### Variable Importance
    if(importance) {
      cat('Refitting models for RF variable importance\n') 
      
      imp_rf_mu_0 <- 
        ranger(outcome_formula, 
               num.threads = 64,
               max.depth = 10,
               num.trees = 500,
               importance = 'impurity_corrected',
               data = df_train %>% filter(surgery == 0),
               seed = 109)
      
      imp_rf_mu_1 <- 
        ranger(outcome_formula, 
               num.threads = 64,
               max.depth = 10,
               num.trees = 500,
               importance = 'impurity_corrected',
               data = df_train %>% filter(surgery == 1),
               seed = 11031)
      
      imp_rf_pi <- 
        ranger(ipw_formula, 
               num.threads = 64,
               max.depth = 2, ### Cap Depth 
               num.trees = 500,
               probability = T,
               importance = 'impurity_corrected',
               data = df_train,
               seed = 81002)
      
      var_importance <- 
        tibble('split_id' = i,
               'covariate' = c(names(imp_rf_mu_0$variable.importance), 
                               names(imp_rf_mu_1$variable.importance),
                               names(imp_rf_pi$variable.importance)),
               'model' = rep(c('mu_0', 'mu_1', 'pi'), c(length(imp_rf_mu_0$variable.importance),
                                                        length(imp_rf_mu_1$variable.importance),
                                                        length(imp_rf_pi$variable.importance))),
               'scaled_importance' = c(imp_rf_mu_0$variable.importance/sum(imp_rf_mu_0$variable.importance),
                                       imp_rf_mu_1$variable.importance/sum(imp_rf_mu_1$variable.importance),
                                       imp_rf_pi$variable.importance/sum(imp_rf_pi$variable.importance))
        )
      
      df_importance <- 
        df_importance %>% 
        bind_rows(var_importance)
      
    }
  }
  
  
  ### Gather Nuisance functions
  df_trials_elig <- 
    df_trials_elig %>% 
    mutate('muA_hat' = ifelse(surgery == 1, mu1_hat, mu0_hat),
           'muA_hat_sub' = ifelse(surgery == 1, mu1_hat_sub, mu0_hat_sub)) %>% 
    ### Weights/stabilization
    group_by(fold_id) %>% 
    mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-pi_hat),
                                     surgery == 1 ~ 1/pi_hat),
           'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat),
                                      surgery == 1 ~ mean(surgery)/(pi_hat)),
           'w_treatment_sub' = case_when(surgery == 0 ~ 1/(1-pi_hat_sub),
                                         surgery == 1 ~ 1/pi_hat_sub)) %>% 
    ungroup() %>%
    group_by(fold_id, sub_fold_id) %>% 
    mutate('sw_treatment_sub' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat_sub),
                                          surgery == 1 ~ mean(surgery)/(pi_hat_sub))) %>% 
    ungroup() %>%
    ### Trim Weights at 99th quantile
    group_by(surgery, fold_id) %>%
    mutate('w_treatment' = winsorize(w_treatment, q = c(0, 0.99)),
           'sw_treatment' = winsorize(sw_treatment, q = c(0, 0.99))) %>%
    ungroup() %>% 
    group_by(surgery, fold_id, sub_fold_id) %>%
    mutate('w_treatment_sub' = winsorize(w_treatment, q = c(0, 0.99)),
           'sw_treatment_sub' = winsorize(sw_treatment, q = c(0, 0.99))) %>%
    ungroup()
  
  ### Trial Specific Chi Hat
  df_trials_elig <- 
    df_trials_elig %>% 
    mutate('chi_IF' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muA_hat),
           'chi_IF_sub' = mu1_hat_sub - mu0_hat_sub + (2 * surgery - 1) * w_treatment_sub * (pct_weight_change - muA_hat_sub))
  
  ### If a subject is eligible in >= 1 trial they have already been assigned a fold 
  ### (e.g. for estimating population eligible fraction) 
  ### If they aren't, we randomize the reamining subjects equally to the folds
  df_index <- 
    df_index %>% 
    left_join(distinct(df_trials_elig, subject_id, fold_id),
              by = 'subject_id')
  
  df_never <- 
    df_index %>% 
    filter(is.na(fold_id)) %>% 
    distinct(subject_id) 
  
  never_split_ids <- 
    split(sample( 1:nrow(df_never) ), ceiling(1:nrow(df_never) / ceiling(nrow(df_never)/n_splits)))
  
  df_never$fold_ix <- 
    map_dfr(1:n_splits, ~{tibble('fold_id' = .x, 'row_ix' = never_split_ids[[.x]])}) %>% 
    arrange(row_ix) %>% 
    pull(fold_id)
  
  df_index <- 
    df_index %>% 
    left_join(df_never, by = 'subject_id') %>% 
    mutate('fold_id' = replace(fold_id, is.na(fold_id), fold_ix[is.na(fold_id)])) %>% 
    select(-fold_ix)
  
  df_IF <- 
    df_trials_elig %>% 
    select(subject_id, trial_id, eligible, surgery, fold_id, sub_fold_id, 
           mu1_hat, mu0_hat, pi_hat, chi_IF,
           mu1_hat_sub, mu0_hat_sub, pi_hat_sub, chi_IF_sub) %>% 
    right_join(df_index, by = c('subject_id', 'trial_id', 'fold_id')) %>% 
    mutate('eligible' = as.numeric(replace(eligible, is.na(eligible), F)),
           'chi_IF' = replace(chi_IF, is.na(chi_IF), 0),
           'chi_IF_sub' = replace(chi_IF_sub, is.na(chi_IF_sub), 0))
  
  
  df_chi <- 
    df_trials_elig %>% 
    group_by(trial_id) %>%
    summarise('chi_DR' = mean(chi_IF),
              'p_treatment' = mean(surgery))
  
  df_info <- 
    df_chi %>% 
    inner_join(df_elig, by = 'trial_id') %>% 
    mutate('outcome' = outcome,
           'procedures' = paste(sort(procedures), collapse = '/'))
  
  if(importance) {
    df_importance <- 
      df_importance %>% 
      group_by(covariate, model) %>% 
      summarise('scaled_importance' = mean(scaled_importance)) %>% 
      ungroup() %>% 
      mutate('outcome' = outcome,
             'procedures' = paste(sort(procedures), collapse = '/'))
  }
  
  return( list('df_info' = df_info,
               'df_importance' = df_importance,
               'df_IF' = df_IF) )
}


### Run analysis for each outcome/procedure(s) combination
df_results <- NULL
df_importance <- NULL
df_seeds <- 
  crossing('outcome' = c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr'),
           'procedures' = c('ALL', 'RYGB', 'SLEEVE')) %>% 
  mutate('seed' = sample(1:100000000, 12))
for(outcome in c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr')) {
  cat('\nAnalysis for Outcome:', outcome, '\n')
  for(procedures in c('ALL', 'RYGB', 'SLEEVE')) {
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
                  ipw_formula = 
                    surgery ~ 
                    baseline_bmi + gender + race + site + baseline_age + 
                    t2dm + insulin + hypertension + hypertension_rx + 
                    dyslipidemia + antilipemic_rx + smoking_status + trial_id,
                  
                  outcome_formula = 
                    pct_weight_change ~
                    baseline_bmi + gender + race + site + baseline_age + 
                    t2dm + insulin + hypertension + hypertension_rx + 
                    dyslipidemia + antilipemic_rx + smoking_status + trial_id,
                  
                  transport_formula = 
                    as.factor(trial_id) ~
                    baseline_bmi + gender + race + site + baseline_age + 
                    t2dm + insulin + hypertension + hypertension_rx + 
                    dyslipidemia + antilipemic_rx + smoking_status,
                  
                  outcome = outcome,
                  n_splits = 2, 
                  n_subsplits = 2, 
                  procedures = procedures,
                  importance = save_importance,
                  trial_range = t_range)
    
    ### Save Results
    df_results <- bind_rows(df_results, pkg$df_info)
    df_importance <-  bind_rows(df_importance, pkg$df_importance)
    
    if(!dir.exists(glue('{data_dir}/tv_effects/IF_contributions/'))) {
      dir.create(glue('{data_dir}/tv_effects/IF_contributions/'))
    }
    
    write_parquet(pkg$df_IF, glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
    
    cat('\n\n---------------------------------\n\n')
  }
}

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
                ipw_formula = 
                  surgery ~ 
                  baseline_bmi + gender + race + site + baseline_age + 
                  t2dm + insulin + hypertension + hypertension_rx + 
                  dyslipidemia + antilipemic_rx + smoking_status + trial_id,
                
                outcome_formula = 
                  pct_weight_change ~
                  baseline_bmi + gender + race + site + baseline_age + 
                  t2dm + insulin + hypertension + hypertension_rx + 
                  dyslipidemia + antilipemic_rx + smoking_status + trial_id,
                
                transport_formula = 
                  as.factor(trial_id) ~
                  baseline_bmi + gender + race + site + baseline_age + 
                  t2dm + insulin + hypertension + hypertension_rx + 
                  dyslipidemia + antilipemic_rx + smoking_status,
                
                outcome = outcome,
                n_splits = 2, 
                n_subsplits = 2, 
                procedures = c('RYGB', 'SLEEVE'),
                importance = save_importance,
                trial_range = 30:84)
  
  ### Save Results
  df_results <- bind_rows(df_results, pkg$df_info)
  df_importance <-  bind_rows(df_importance, pkg$df_importance)
  
  if(!dir.exists(glue('{data_dir}/tv_effects/IF_contributions/'))) {
    dir.create(glue('{data_dir}/tv_effects/IF_contributions/'))
  }
  
  write_parquet(pkg$df_IF, glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome}_{gsub("/", "-", pkg$df_info$procedures[1])}.parquet'))
  
  cat('\n\n---------------------------------\n\n')
}

### Save Results
write_csv(df_results, glue('{data_dir}/tv_effects/weight_DR_results.csv'))
if(save_importance) {
  write_csv(df_importance, glue('{data_dir}/tv_effects/weight_DR_importance.csv'))
}

