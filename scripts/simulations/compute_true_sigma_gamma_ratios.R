library(tidyverse)
source('scripts/simulations/generate_data.R') 
source('scripts/simulations/sim_helpers.R') 

files <- dir('sim_inputs')
files <- files[grepl('params', files)]
set.seed(1791317)
df_summary <- NULL
n_rep <- 10
for(rep in 1:n_rep) {
  for(i in 1:length(files)) {
    cat('Rep', rep,'of', n_rep, '; Combo', i, 'of', length(files), '\n')
    params <- read_rds(glue('sim_inputs/params_{i}.rds'))
    
    params$sigma_rss <- 1e-10
    params$n_subjects <- 1000
    if(params$shift != 'No Covariate Shift' & grepl('Effect Modification', params$trt_change)) {
      params$n_subjects <- 50000 
    }
    
    df_trials <- generate_data(params)
    
    
    df_trials_elig <- 
      df_trials %>% 
      filter(eligible)
    
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
      ) 
    
    
    ### Get the "True" Values of mu, pi
    baseline_trials <- 
      baseline_trials %>% 
      mutate('mu1_hat' = compute_model(baseline_trials %>% mutate('surgery' = 1), params$outcome_model),
             'mu0_hat' = compute_model(baseline_trials %>% mutate('surgery' = 0), params$outcome_model),
             'pi_hat' = expit(compute_model(baseline_trials, params$treatment_model))) %>% 
      group_by(subject_id, baseline_trial_id) %>% 
      mutate('mu1j_hat' = mu1_hat[trial_id == baseline_trial_id],
             'mu0j_hat' = mu0_hat[trial_id == baseline_trial_id]) %>% 
      ungroup() %>% 
      mutate('muA_hat' = ifelse(surgery == 1, mu1_hat, mu0_hat),
             'muAj_hat' = ifelse(surgery == 1, mu1j_hat, mu0j_hat)) %>% 
      mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-pi_hat),
                                       surgery == 1 ~ 1/pi_hat),
             'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-pi_hat),
                                        surgery == 1 ~ mean(surgery)/(pi_hat)))
    
    standardize_IF <- 
      baseline_trials %>% 
      mutate('chi_IF' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muA_hat),
             'chi_IFj' = mu1_hat - mu0_hat + (2 * surgery - 1) * w_treatment * (pct_weight_change - muAj_hat)) %>% 
      mutate('eligible' = as.numeric(replace(eligible, is.na(eligible), F)),
             'chi_IF' = replace(chi_IF, is.na(chi_IF), 0),
             'chi_IFj' = replace(chi_IFj, is.na(chi_IFj), 0))
    
    ### Standarization Matrix
    standardization_matrix <- 
      standardize_IF %>% 
      group_by(baseline_trial_id, trial_id) %>%
      summarise('chi_DR' = mean(chi_IF),
                'chi_DRj' = mean(chi_IFj),
                'chi_gformula' = mean(mu1_hat - mu0_hat)) %>% 
      ungroup()
    
    df_tmp <- 
      sm_summary_metrics(standardization_matrix) %>% 
      mutate('trt_change' = params$trt_change,
             'shift' = params$shift)
    
    df_summary <- 
      df_summary %>% 
      bind_rows(df_tmp)
  }
}

df_summary <- 
  df_summary %>% 
  group_by(trial_id, estimator, shift, trt_change) %>% 
  summarise('gamma2_m' = mean(gamma2_m),
            'sigma2_m' = mean(sigma2_m),
            'sigma_ratio' = mean(sigma_ratio),
            'sigma2_ratio' = mean(sigma2_ratio)) %>% 
  ungroup()

write_csv(df_summary, 'sim_inputs/true_sigma_ratios.csv') 
