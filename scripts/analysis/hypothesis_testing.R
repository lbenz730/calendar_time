library(tidyverse)
library(glue)
library(arrow)
library(data.table)
library(matrixStats)

source('scripts/simulations/sim_helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Standardization Analysis We already Computed
df_transport <- 
  map_dfr(1:84, ~{
    read_csv(glue('{data_dir}/tv_effects/weight_DR_results_baseline_trial_{.x}.csv')) %>% 
      mutate('baseline_trial' = .x)
  }) 

load_IF <- function(baseline_trial_id, outcome, procedures) {
  ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
  data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
  
  file_path <- glue('{data_dir}/tv_effects/IF_contributions/baseline_trial_{baseline_trial_id}/subject_IF_cross_{outcome}_{gsub("/", "-", procedures)}.parquet')
  df_IF <- read_parquet(file_path)
  
  if(procedures %in% c('SLEEVE', 'RYGB/SLEEVE')) {
   df_IF <- 
     df_IF %>% 
     filter(trial_id >= 30)
  }
  
  df_IF <-
    df_IF %>% 
    mutate('subject_IF' = IF_cross) %>% 
    mutate('outcome' = outcome, 
           'procedures' = procedures) %>% 
    select(outcome, procedures, subject_id, trial_id, subject_IF) %>% 
    pivot_wider(names_from = 'trial_id',
                values_from = 'subject_IF', 
                names_prefix = 'm_') %>% 
    set_names(c('outcome', 'procedures', 'subject_id'), paste0(names(.)[-c(1:3)], '_j_', baseline_trial_id)) %>% 
    arrange(subject_id) 
  
  return(df_IF)
}


hypothesis_test <- function(outcome_, procedures_, n_boot = 10000) {

  ### Standardization Matrix for Current Analysis
  standardization_df <- 
    df_transport %>% 
    filter(outcome == outcome_, 
           procedures == procedures_) 
  
  M <- n_distinct(standardization_df$trial_id)
  standardization_matrix <- matrix(standardization_df$chi_cross, nrow = M, ncol = M, byrow = T)
  
  chi_vec <- diag(standardization_matrix)
  row_null <- outer(chi_vec, rep(1, M))
  col_null <- t(outer(chi_vec, rep(1, M)))
  sigma2 <- rowSums( (standardization_matrix - row_null )^2 )/(M-1) 
  gamma2 <- colSums( (standardization_matrix -  col_null)^2 )/(M-1) 
  test_stat <- mean(sigma2/(sigma2 + gamma2))
  
  ### Load in All IF Contributions to Estimate Covariance Matrix
  cat('Loading IF Contributions\n')
  df_IF <- 
    map(sort(unique(standardization_df$trial_id)), ~{
      cat('Loading', .x, '\n')
      load_IF(baseline_trial_id = .x,
              outcome = outcome_,
              procedures = procedures_)
    })
  
  ### Convert to DT for efficient combination
  cat('Merge IF Contributions\n')
  df_IF <- lapply(df_IF, as.data.table)
  df_IF <- lapply(df_IF, setkey, outcome, procedures, subject_id)
  
  ### Iteratively merge with data.table (fast + low memory overhead)
  df_IF <- Reduce(function(x, y) x[y, nomatch = 0], df_IF)
  
  ### Estimate Covariance Matrix
  cat('Estimate Covariance Matrix (Matrix Conversion)\n')
  X_IF <- as.matrix(df_IF[,-c(1:3)])
  # X_IF <- scale(X_IF, center = T, scale = F)   # center columns
  cat('Estimate Covariance Matrix (Center Cols)\n')
  X_IF <- sweep(X_IF, 2, colMeans(X_IF), "-") # center columns
  cat('Estimate Covariance Matrix (Cross Product)\n')
  Sigma_IF <- ( crossprod(X_IF) / (nrow(X_IF) - 1) ) * 1/nrow(X_IF)
  
  ### Parametric Bootstrap
  cat('Parametric Bootstrap\n')
  resamples <- MASS::mvrnorm(n = n_boot, mu = as.vector(t(standardization_matrix)), Sigma = Sigma_IF) 
  sim_stats <- map_dbl(1:n_boot, ~fast_ratio(vec = resamples[.x,], M, threshold = 0))
  q <- quantile(sim_stats, c(0.025, 0.975))
  
  df_test <- 
    tibble('outcome' = outcome_, 
           'procedures' = procedures_,
           'test_stat' = test_stat,
           'sd_stat' = sd(sim_stats),
           'qlow' = q[1],
           'qhigh' = q[2])
  
  return(df_test)
  
}

### Conduct Tests
set.seed(617)
df_hypothesis_all <- NULL
df_seeds <- 
  crossing('outcome' = c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr'),
           'procedures' = c('AGB/RYGB/SLEEVE', 'RYGB', 'SLEEVE', 'RYGB/SLEEVE')) %>% 
  mutate('seed' = sample(1:100000000, 16))

for(setting in 1:16) {
  cat('Setting:', setting, '\n')
  set.seed(df_seeds$seed[setting])
  
  df_hypothesis <- 
    hypothesis_test(outcome_ = df_seeds$outcome[setting], 
                    procedures_ = df_seeds$procedures[setting],
                    n_boot = 10000)
  
  write_csv(df_hypothesis, glue('{data_dir}/tv_effects/EIF_weight_hypothesis_tests_{setting}.csv'))
  
  gc()
}