### Wrapper to run one simulation iteration for all parameters
### So that we can parallelize across iterations rather than settings

library(tidyverse)
library(glue)
library(arrow)
library(splines)
library(ranger)
library(nleqslv)
library(SuperLearner)

source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R') 
source('scripts/simulations//sim_helpers.R') 
source('scripts/simulations/generate_data.R')

### Simulation Settings
args <- commandArgs(trailingOnly = T)
sim_id <- as.numeric(args[1])
n_subj <- as.numeric(args[2]) ### number of subjects to use
n_threads <- ifelse(n_subj < 10000, 16, 32)
n_sims <- 1000
n_params <- 18
equiv_threshold <- 0.005

### Iteration specific seed
set.seed(123)
seeds <- sample(1:1000000000, n_sims)
set.seed(seeds[sim_id])

### Different Seed for one iteration of small simulation 
### to avoid splits w/ no factor levels 
if(sim_id %in% c(424, 653, 835) & n_subj == 100) {
  set.seed(seeds[sim_id] + 1)
}

### Chi DR EIF Function
### Also saves standardization matrix
fit_EIF_chi <- function(df_trials, mu_covariates, pi_covariates, xi_covariates, n_splits, n_subsplits) {
  
  cat('Computing Empirical Eligibility\n')
  ### Empirical Eligibility
  df_trials_elig <- 
    df_trials %>% 
    filter(eligible)
  
  df_elig <- 
    df_trials %>% 
    mutate('eligible' = ifelse(is.na(eligible), F, eligible)) %>% 
    group_by(trial_id) %>% 
    summarise('p_elig' = mean(eligible)) 
  
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
  
  
  ### Storage for Standardization Analysis
  baseline_trials <- 
    df_trials_elig %>% 
    rename('baseline_trial_id' = trial_id) %>% 
    inner_join(
      crossing('subject_id' = unique(df_trials_elig$subject_id),
               'baseline_trial_id' = unique(df_trials_elig$trial_id),
               'trial_id' = unique(df_trials_elig$trial_id)),
      by = c('subject_id', 'baseline_trial_id')
    ) %>% 
    ### Update covariates that depend on treatment time
    select(-contains('spline')) %>% 
    inner_join(
      as_tibble(ns(1:max(df_trials_elig$trial_id), df = 3)) %>% 
        mutate_at(vars(everything()), as.numeric) %>% 
        set_names(paste0('spline_', 1:ncol(.))) %>% 
        mutate('trial_id' = 1:max(df_trials_elig$trial_id)),
      by = 'trial_id'
    ) %>% 
    group_by(subject_id, trial_id) %>% 
    mutate('surgery_m' = surgery[trial_id == baseline_trial_id],
           'pct_weight_change_m' = pct_weight_change[trial_id == baseline_trial_id]) %>% 
    ungroup()
  
  baseline_trials <- 
    baseline_trials %>% 
    mutate('fold_id' = NA_real_,
           'row_ix' = 1:nrow(.),
           'mu0_hat' = NA_real_,
           'mu1_hat' = NA_real_,
           'pi_hat' = NA_real_,
           'pim_hat' = NA_real_,
           'mu0j_hat' = NA_real_,
           'mu1j_hat' = NA_real_,
           'mu0m_hat' = NA_real_,
           'mu1m_hat' = NA_real_,
           'xi_m' = NA_real_,
           'xi_j' = NA_real_)
  
  ### Set Up Ranger RF Learners for SuperLearner
  rf_lnr_mu <- 
    create.Learner('SL.ranger', tune = list('max.depth' = 10,
                                            'num.trees' = 500,
                                            'num.threads' = n_threads))
  
  rf_lnr_pi <- 
    create.Learner('SL.ranger', tune = list('max.depth' = 4,
                                            'num.trees' = 500,
                                            'num.threads' = n_threads))
  
  ### Set up GAM learner
  gam_lnr <- 
    create.Learner('SL.gam', tune = list('deg.gam' = 4))
  
  
  for(i in 1:n_splits) {
    test_ids <- subj_ids[ split_ids[[i]] ]
    train_ids <- subj_ids[ setdiff(1:length(subj_ids), split_ids[[i]]) ]
    
    df_train <- 
      df_trials_elig %>% 
      filter(subject_id %in% train_ids)
    
    df_test <- 
      df_trials_elig %>% 
      filter(subject_id %in% test_ids)
    
    ### Hold out predictions for baseline trial population
    df_test_baseline <- 
      baseline_trials %>% 
      filter(subject_id %in% test_ids)
    
    ### Outer Level Fitting 
    cat('Fitting Outcome Model on split [', i, '/', n_splits, ']\n', sep = '')
    ### Fit Outcome Model
    
    Y_train_0 <- df_train$pct_weight_change[df_train$surgery == 0]
    Y_train_1 <- df_train$pct_weight_change[df_train$surgery == 1]
    
    X_train_0 <- 
      df_train %>% 
      filter(surgery == 0) %>% 
      select(all_of(mu_covariates)) 
    
    X_train_1 <- 
      df_train %>% 
      filter(surgery == 1) %>% 
      select(all_of(mu_covariates)) 
    
    X_test_mu <- 
      df_test %>% 
      select(all_of(mu_covariates))
    
    ### Fit Super Learner for Outcome Model
    sl_mu_0 <- 
      SuperLearner(Y = Y_train_0,
                   X = one_hot_encode(X_train_0), 
                   cvControl = list(V = 3),
                   SL.library = c(rf_lnr_mu$names, gam_lnr$names, 'SL.speedglm'))
    
    sl_mu_1 <- 
      SuperLearner(Y = Y_train_1,
                   X = one_hot_encode(X_train_1), 
                   cvControl = list(V = 3),
                   SL.library = c(rf_lnr_mu$names, gam_lnr$names, 'SL.speedglm'))
    
    cat('Fitting PS Model on split [', i, '/', n_splits, ']\n', sep = '')
    Y_train <- df_train$surgery
    X_train <- 
      df_train %>%   
      select(all_of(pi_covariates))
    
    X_test_pi <- 
      df_test %>% 
      select(all_of(pi_covariates))
    
    sl_pi <- 
      SuperLearner(Y = Y_train,
                   X = one_hot_encode(X_train), 
                   family = 'binomial',
                   cvControl = list(V = 3),
                   SL.library = c(rf_lnr_pi$names, gam_lnr$names, 'SL.glm'))
    
    
    cat('Fitting Transport Model on split [', i, '/', n_splits, ']\n', sep = '')
    Y_train <- as.factor(df_train$trial_id)
    X_train <- 
      df_train %>% 
      select(trial_id, all_of(xi_covariates))
    
    rf_xi <- 
      ranger(trial_id ~ ., 
             num.threads = n_threads,
             max.depth = 10,
             num.trees = 500,
             data = X_train, 
             probability = T)
    
    
    ### Save Predictions and Models
    cat('Predicting mu0, mu1, pi models on test set [', i, '/', n_splits, ']\n', sep = '')
    df_test$mu0_hat <- SL.predict(sl_mu_0, clean_test_df(X_test_mu, X_train_0))
    df_test$mu1_hat <- SL.predict(sl_mu_1, clean_test_df(X_test_mu, X_train_1))
    df_test$pi_hat <- SL.predict(sl_pi, clean_test_df(X_test_pi, X_train), binomial = T)
    
    ### Transport Matrix for Density Ratios
    cat('Predicting Transport Matrix for Density Ratios [', i, '/', n_splits, ']\n', sep = '')
    transport_matrix <- predict(rf_xi, data = df_test)$predictions
    colnames(transport_matrix) <- paste0('xi_', sort(unique(df_trials$trial_id)))
    transport_matrix <- 
      df_test %>% 
      select(subject_id, trial_id) %>% 
      bind_cols(as_tibble(transport_matrix))
    
    transport_df <- 
      transport_matrix %>% 
      group_by(subject_id, trial_id) %>% 
      pivot_longer(cols = contains('xi_'),
                   names_to = 'baseline_trial_id', 
                   values_to = 'xi_j',
                   names_prefix = 'xi_') %>% 
      mutate('baseline_trial_id' = as.numeric(baseline_trial_id)) %>% 
      group_by(subject_id, trial_id) %>% 
      mutate('xi_m' = xi_j[trial_id == baseline_trial_id]) %>% 
      ungroup()
    
    df_trials_elig$mu0_hat[df_test$row_ix] <- df_test$mu0_hat 
    df_trials_elig$mu1_hat[df_test$row_ix] <- df_test$mu1_hat 
    df_trials_elig$pi_hat[df_test$row_ix] <- df_test$pi_hat 
    df_trials_elig$fold_id[df_test$row_ix] <- i ### Save the index of the holdout fold(s)
    
    X_test_baseline_mu <- 
      df_test_baseline %>% 
      select(all_of(mu_covariates), baseline_trial_id)
    
    X_test_baseline_pi <- 
      df_test_baseline %>% 
      select(all_of(pi_covariates), baseline_trial_id)
    
    df_test_baseline$mu0_hat <- SL.predict(sl_mu_0, clean_test_df(X_test_baseline_mu %>% select(-baseline_trial_id), X_train_0))
    df_test_baseline$mu1_hat <- SL.predict(sl_mu_1, clean_test_df(X_test_baseline_mu %>% select(-baseline_trial_id), X_train_1))
    df_test_baseline$mu0j_hat <- SL.predict(sl_mu_0, clean_test_df(X_test_baseline_mu %>% mutate('trial_id' = baseline_trial_id) %>% select(-baseline_trial_id), X_train_0))
    df_test_baseline$mu1j_hat <- SL.predict(sl_mu_1, clean_test_df(X_test_baseline_mu %>% mutate('trial_id' = baseline_trial_id) %>% select(-baseline_trial_id), X_train_1))
    df_test_baseline$pi_hat <- SL.predict(sl_pi, clean_test_df(X_test_baseline_pi %>% select(-baseline_trial_id), X_train), binomial = T)
    
    df_test_baseline <- 
      df_test_baseline %>% 
      select(-contains('xi_')) %>% 
      inner_join(transport_df, by = c('subject_id', 'trial_id', 'baseline_trial_id')) 
    
    df_test_baseline <- 
      df_test_baseline %>% 
      group_by(subject_id, trial_id) %>% 
      mutate('mu1m_hat' = mu1j_hat[trial_id == baseline_trial_id],
             'mu0m_hat' = mu0j_hat[trial_id == baseline_trial_id],
             'pim_hat' = pi_hat[trial_id == baseline_trial_id])
    
    baseline_trials$mu0_hat[df_test_baseline$row_ix] <- df_test_baseline$mu0_hat 
    baseline_trials$mu1_hat[df_test_baseline$row_ix] <- df_test_baseline$mu1_hat 
    baseline_trials$mu0j_hat[df_test_baseline$row_ix] <- df_test_baseline$mu0j_hat 
    baseline_trials$mu1j_hat[df_test_baseline$row_ix] <- df_test_baseline$mu1j_hat 
    baseline_trials$mu0m_hat[df_test_baseline$row_ix] <- df_test_baseline$mu0m_hat 
    baseline_trials$mu1m_hat[df_test_baseline$row_ix] <- df_test_baseline$mu1m_hat 
    baseline_trials$pi_hat[df_test_baseline$row_ix] <- df_test_baseline$pi_hat 
    baseline_trials$pim_hat[df_test_baseline$row_ix] <- df_test_baseline$pim_hat 
    baseline_trials$xi_m[df_test_baseline$row_ix] <- df_test_baseline$xi_m
    baseline_trials$xi_j[df_test_baseline$row_ix] <- df_test_baseline$xi_j
    baseline_trials$fold_id[df_test_baseline$row_ix] <- i
    
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
      ### Fit Outcome Model
      
      Y_train_sub_0 <- df_train_sub$pct_weight_change[df_train_sub$surgery == 0]
      Y_train_sub_1 <- df_train_sub$pct_weight_change[df_train_sub$surgery == 1]
      
      X_train_sub_0 <- 
        df_train_sub %>% 
        filter(surgery == 0) %>% 
        select(all_of(mu_covariates))
      
      X_train_sub_1 <- 
        df_train_sub %>% 
        filter(surgery == 1) %>% 
        select(all_of(mu_covariates))
      
      X_test_sub_mu <- 
        df_test_sub %>% 
        select(all_of(mu_covariates))
      
      ### Fit Super Learner for Outcome Model
      sl_mu_0_sub <- 
        SuperLearner(Y = Y_train_sub_0,
                     X = one_hot_encode(X_train_sub_0), 
                     cvControl = list(V = 3),
                     SL.library = c(rf_lnr_mu$names, gam_lnr$names, 'SL.speedglm'))
      sl_mu_1_sub <- 
        SuperLearner(Y = Y_train_sub_1,
                     X = one_hot_encode(X_train_sub_1), 
                     cvControl = list(V = 3),
                     SL.library = c(rf_lnr_mu$names, gam_lnr$names, 'SL.speedglm'))
      
      ### Fit Propensity Model
      cat('Fitting PS Model on subsplit [', j, '/', n_subsplits, '] for training fold ', i, '\n', sep = '')
      Y_train_sub <- df_train_sub$surgery
      X_train_sub <- 
        df_train_sub %>%   
        select(all_of(pi_covariates))
      
      X_test_sub_pi <- 
        df_test_sub %>% 
        select(all_of(pi_covariates))
      
      sl_pi_sub <- 
        SuperLearner(Y = Y_train_sub,
                     X = one_hot_encode(X_train_sub), 
                     family = 'binomial',
                     cvControl = list(V = 3),
                     SL.library = c(rf_lnr_pi$names, 'SL.gam', 'SL.glm'))
      
      ### Save Predictions
      cat('Predicting mu0, mu1, pi models on test subsplit [', j, '/', n_subsplits, '] for training fold ', i, '\n', sep = '')
      df_test_sub$mu0_hat <- SL.predict(sl_mu_0_sub, clean_test_df(X_test_sub_mu, X_train_sub_0))
      df_test_sub$mu1_hat <- SL.predict(sl_mu_1_sub, clean_test_df(X_test_sub_mu, X_train_sub_1))
      df_test_sub$pi_hat <- SL.predict(sl_pi_sub, clean_test_df(X_test_sub_pi, X_train_sub), binomial = T)
      
      df_trials_elig$mu0_hat_sub[df_test_sub$row_ix] <- df_test_sub$mu0_hat 
      df_trials_elig$mu1_hat_sub[df_test_sub$row_ix] <- df_test_sub$mu1_hat 
      df_trials_elig$pi_hat_sub[df_test_sub$row_ix] <- df_test_sub$pi_hat 
      df_trials_elig$sub_fold_id[df_test_sub$row_ix] <- j ### Save the index of the holdout fold(s)
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
  
  weight_max <- 
    df_trials_elig %>% 
    group_by(surgery, fold_id) %>%
    summarise('w_treatment_bound' = max(w_treatment),
              'sw_treatment_bound' = max(sw_treatment)) %>%
    ungroup()
  
  weight_max_m <- 
    df_trials_elig %>% 
    group_by('surgery_m' = surgery, fold_id) %>%
    summarise('w_treatment_bound_m' = max(w_treatment),
              'sw_treatment_bound_m' = max(sw_treatment)) %>%
    ungroup()
  
  ### Trial Specific Chi Hat
  df_trials_elig <- 
    df_trials_elig %>% 
    mutate('chi_IF' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muA_hat),
           'chi_IF_sub' = mu1_hat_sub - mu0_hat_sub + (2 * surgery - 1) * w_treatment_sub * (pct_weight_change - muA_hat_sub))
  
  df_IF <- 
    df_trials_elig %>% 
    select(subject_id, trial_id, eligible, surgery, fold_id, sub_fold_id, 
           mu1_hat, mu0_hat, pi_hat, chi_IF,
           mu1_hat_sub, mu0_hat_sub, pi_hat_sub, chi_IF_sub) %>% 
    mutate('eligible' = as.numeric(replace(eligible, is.na(eligible), F)),
           'chi_IF' = replace(chi_IF, is.na(chi_IF), 0),
           'chi_IF_sub' = replace(chi_IF_sub, is.na(chi_IF_sub), 0))
  
  
  ### Standardization Analysis
  baseline_trials <- 
    baseline_trials %>% 
    mutate('muA_hat' = ifelse(surgery == 1, mu1_hat, mu0_hat),
           'muAj_hat' = ifelse(surgery == 1, mu1j_hat, mu0j_hat),
           'muAm_hat' = ifelse(surgery_m == 1, mu1m_hat, mu0m_hat)) %>% 
    ### Weights/stabilization
    group_by(surgery, fold_id) %>%
    mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-pi_hat),
                                     surgery == 1 ~ 1/pi_hat),
           'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat),
                                      surgery == 1 ~ mean(surgery)/(pi_hat)),
           'w_treatment_m' = case_when(surgery_m == 0 ~ 1/(1-pim_hat),
                                       surgery_m == 1 ~ 1/pim_hat),
           'sw_treatment_m' = case_when(surgery_m == 0 ~ mean(1-surgery_m)/(1-pim_hat),
                                        surgery_m == 1 ~ mean(surgery_m)/(pim_hat))) %>% 
    inner_join(weight_max, by = c('surgery', 'fold_id')) %>% 
    inner_join(weight_max_m, by = c('surgery_m', 'fold_id')) %>% ### here edit
    ### Trim Weights at 99th quantile (overall)
    mutate('w_treatment' = pmin(w_treatment, w_treatment_bound),
           'sw_treatment' = pmin(sw_treatment, sw_treatment_bound),
           'w_treatment_m' = pmin(w_treatment_m, w_treatment_bound_m),
           'sw_treatment_m' = pmin(sw_treatment_m, sw_treatment_bound_m)) %>% 
  ungroup()

standardize_IF <- 
  baseline_trials %>% 
  mutate('chi_IF' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muA_hat),
         'chi_IFj' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muAj_hat),
         'chi_cross_IF' = mu1_hat - mu0_hat + (2 * surgery_m - 1) * w_treatment_m *  xi_j/xi_m * (pct_weight_change_m - muAm_hat)) %>% 
  select(subject_id, baseline_trial_id, trial_id, fold_id, eligible, 
         surgery, surgery_m, pct_weight_change, pct_weight_change_m,
         mu1_hat, mu0_hat, muA_hat, mu1j_hat, mu0j_hat, muAj_hat, mu1m_hat, mu0m_hat, muAm_hat, 
         pi_hat, pim_hat, w_treatment, w_treatment_m, xi_m, xi_j,
         chi_IF, chi_IFj, chi_cross_IF) %>% 
  mutate('eligible' = as.numeric(replace(eligible, is.na(eligible), F)),
         'chi_IF' = replace(chi_IF, is.na(chi_IF), 0),
         'chi_IFj' = replace(chi_IFj, is.na(chi_IFj), 0))

standardization_matrix <- 
  standardize_IF %>% 
  group_by(baseline_trial_id, trial_id) %>%
  summarise('chi_DR' = mean(chi_IF),
            'chi_DRj' = mean(chi_IFj),
            'chi_cross' = mean(chi_cross_IF),
            'chi_gformula' = mean(mu1_hat - mu0_hat),
            'p_treatment' = mean(surgery)) %>% 
  ungroup()

return( list('df_IF' = df_IF, 
             'standardize_IF' = standardize_IF,
             'standardization_matrix' = standardization_matrix,
             'baseline_trials' = baseline_trials) )
}


### Loop over parameter settings to do a single simulation iteration for each
df_metrics_all <- NULL
df_loss_all <- NULL
df_hypothesis_all <- NULL
for(setting in 1:n_params) {
  cat('SIMULATION SETTING:', setting, 'OF', n_params, '\n')
  ### Load parameters 
  params <- read_rds(glue('sim_inputs/params_{setting}.rds'))
  params$n_subjects <- n_subj
  params$equiv_threshold <- equiv_threshold
  
  ### Generate Simulated Dataset
  df_trials <- generate_data(params)
  
  ### Run analysis to get EIF Contributions and Standardization Matrix
  pkg <- 
    fit_EIF_chi(df_trials = df_trials, 
                mu_covariates =   
                  c('baseline_bmi', 'gender', 'race', 'site', 'baseline_age',
                    't2dm', 'insulin', 'hypertension', 'hypertension_rx',
                    'dyslipidemia', 'antilipemic_rx', 'smoking_status', 'trial_id'),
                pi_covariates =   
                  c('baseline_bmi', 'gender', 'race', 'site', 'baseline_age',
                    't2dm', 'insulin', 'hypertension', 'hypertension_rx',
                    'dyslipidemia', 'antilipemic_rx', 'smoking_status', 'trial_id'),
                xi_covariates =   
                  c('baseline_bmi', 'gender', 'race', 'site', 'baseline_age',
                    't2dm', 'insulin', 'hypertension', 'hypertension_rx',
                    'dyslipidemia', 'antilipemic_rx', 'smoking_status'),
                n_splits = 2, 
                n_subsplits = 2)
  
  ### Projection Step
  df_loss <- 
    EIF_projection(pkg$df_IF) %>% 
    mutate('parameter_id' = setting,
           'sim_id' = sim_id,
           'n_trials' = params$n_trials,
           'n_subjects' = params$n_subjects,
           'shift' = params$shift,
           'trt_change' = params$trt_change)
  
  ### Compute Standardization Matrix summary statistics 
  df_metrics <- 
    sm_summary_metrics(pkg$standardization_matrix, threshold = params$equiv_threshold) %>% 
    mutate('parameter_id' = setting,
           'sim_id' = sim_id,
           'n_trials' = params$n_trials,
           'n_subjects' = params$n_subjects,
           'equiv_threshold' = params$equiv_threshold,
           'shift' = params$shift,
           'trt_change' = params$trt_change)
  
  ### Hypothesis Testing
  df_hypothesis <- 
    boot_ratio_ci(standardize_IF = pkg$standardize_IF, 
                  standardization_df = pkg$standardization_matrix, 
                  threshold = params$equiv_threshold,
                  M = params$n_trials,
                  n_boot = 10000) %>% 
    mutate('parameter_id' = setting,
           'sim_id' = sim_id,
           'n_trials' = params$n_trials,
           'n_subjects' = params$n_subjects,
           'shift' = params$shift,
           'trt_change' = params$trt_change)
  
  df_metrics_all <- 
    df_metrics_all %>% 
    bind_rows(df_metrics)
  
  df_hypothesis_all <- 
    df_hypothesis_all %>% 
    bind_rows(df_hypothesis) 
  
  df_loss_all <- 
    df_loss_all %>% 
    bind_rows(df_loss) 
  
  cat('---------------------------------\n')
}

### Write Out Iteration Results 
write_parquet(df_metrics_all, glue('sim_outputs/n{n_subj}/standardization_metrics/sim_{sim_id}.parquet'))
write_parquet(df_loss_all, glue('sim_outputs/n{n_subj}/loss_fx/sim_{sim_id}.parquet'))
write_parquet(df_hypothesis_all, glue('sim_outputs/n{n_subj}/hypothesis_testing/sim_{sim_id}.parquet'))
