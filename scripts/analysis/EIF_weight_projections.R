library(tidyverse)
library(glue)
library(arrow)
library(splines)
library(mgcv)
library(nleqslv)
library(ranger)
source('scripts/helpers.R')
source('scripts/analysis/EIF_helpers.R')

set.seed(04071318)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Model Selection of Candidate Projections
df_risk <- NULL
df_mse <- NULL
for(outcome_ in c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr')) {
  for(procedures_ in c('AGB/RYGB/SLEEVE', 'RYGB', 'SLEEVE', 'RYGB/SLEEVE')) {
    cat(outcome_, procedures_, '\n')
    
    if(procedures_ == 'AGB/RYGB/SLEEVE') {
      t_range <- 1:84 
      sk2 <- c(24, 48)
      sk3 <- c(12, 36, 60)
    } else if(procedures_ == 'RYGB') {
      t_range <- 1:84 
      sk2 <- c(24, 48)
      sk3 <- c(12, 36, 60)
    } else if(procedures_ == 'SLEEVE'){
      t_range <- 30:84
      sk2 <- c(48, 72)
      sk3 <- c(40, 56, 72)
    }  else if(procedures_ == 'RYGB/SLEEVE'){
      t_range <- 30:84
      sk2 <- c(48, 72)
      sk3 <- c(40, 56, 72)
    }
    
    ### Subject Specific IF Contributions
    df_IF <- 
      read_parquet(glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome_}_{gsub("/", "-", procedures_)}.parquet')) %>% 
      filter(trial_id %in% t_range)
    
    
    df_chi <-
      df_IF %>%
      group_by(trial_id, fold_id) %>%
      summarise('chi_DR' = mean(chi_IF[eligible == 1]),
                'chi_DR_for_proj' = mean(chi_IF_sub[eligible == 1]),
                'p_treatment' = mean(surgery[eligible == 1]),
                'p_elig' = mean(eligible),
                'n' = n(),
                'n_elig' = sum(eligible)) %>% 
      group_by(trial_id) %>% 
      mutate('p_elig_for_proj' = (sum(n_elig) - n_elig)/(sum(n) - n)) %>% 
      ungroup() %>% 
      select(-n, -n_elig)
    
    ### w(m) options
    ### unifom w(m) = 1
    ### inv_elig = 1/P(E[m] = 1) (eg. to cancel out P(E[m] = 1) in estimating eqn
    for(fold in 1:max(df_chi$fold_id)) {
      df_fold <- 
        df_chi %>% 
        filter(fold_id == fold)
      
      for(weights in c('uniform')) {
        if(weights == 'uniform') {
          w <- rep(1, length(t_range))
        } else if(weights == 'inv_elig') {
          w <- 1/df_info$p_elig 
        }
        
        ### Constant projection
        beta_constant <- 
          nleqslv(x = c(0), 
                  fn = estimating_eqn, 
                  trial_range = t_range,
                  weights = w,
                  form = 'constant', 
                  ate_vec = df_fold$chi_DR_for_proj,
                  p_elig = df_fold$p_elig_for_proj)
        
        ### Linear Projection
        beta_linear <- 
          nleqslv(x = c(0,0), 
                  fn = estimating_eqn, 
                  trial_range = t_range,
                  weights = w,
                  form = 'polynomial', 
                  ate_vec = df_fold$chi_DR_for_proj,
                  p_elig = df_fold$p_elig_for_proj,
                  poly_k = 1)
        
        ### Cubic Projection
        cubic_start <- 
          lm(chi_DR_for_proj ~ trial_id + I(trial_id^2) + I(trial_id^3), 
             data = df_fold,
             weight = p_elig_for_proj)
        
        beta_cubic <- 
          nleqslv(x = coef(cubic_start), 
                  fn = estimating_eqn, 
                  trial_range = t_range,
                  weights = w,
                  form = 'polynomial', 
                  ate_vec = df_fold$chi_DR_for_proj,
                  p_elig = df_fold$p_elig_for_proj,
                  poly_k = 3)
        
        beta_spline_2 <- 
          nleqslv(x = c(0,0,0,0), 
                  fn = estimating_eqn, 
                  trial_range = t_range,
                  weights = w,
                  form = 'spline', 
                  ate_vec = df_fold$chi_DR_for_proj,
                  p_elig = df_fold$p_elig_for_proj,
                  spline_knots = sk2)
        
        beta_spline_3 <- 
          nleqslv(x = c(0,0,0,0,0), 
                  fn = estimating_eqn, 
                  trial_range = t_range,
                  weights = w,
                  form = 'spline', 
                  ate_vec = df_fold$chi_DR_for_proj,
                  p_elig = df_fold$p_elig_for_proj,
                  spline_knots = sk3)
        
        
        ### Design Matricies
        basis2 <- ns(x = t_range, knots = sk2)
        basis3 <- ns(x = t_range, knots = sk3)
        X_linear <- cbind(1, t_range)
        X_cubic <- cbind(1, t_range, t_range^2, t_range^3)
        X_s2 <- cbind(1, basis2)
        X_s3 <- cbind(1, basis3)
        
        ### Projection Variances
        var_beta_constant <-
          var_beta(df_IF = df_IF,
                   beta = beta_constant$x,
                   trial_range = t_range,
                   weights = w,
                   form = 'constant',
                   ate_vec =  df_fold$chi_DR_for_proj,
                   p_elig = df_fold$p_elig_for_proj)
        
        var_beta_linear <-
          var_beta(df_IF = df_IF,
                   beta = beta_linear$x,
                   trial_range = t_range,
                   weights = w,
                   form = 'polynomial',
                   ate_vec =  df_fold$chi_DR_for_proj,
                   p_elig = df_fold$p_elig_for_proj,
                   poly_k = 1)
        
        
        var_beta_cubic <-
          var_beta(df_IF = df_IF,
                   beta = beta_cubic$x,
                   trial_range = t_range,
                   weights = w,
                   form = 'polynomial',
                   ate_vec =  df_fold$chi_DR_for_proj,
                   p_elig = df_fold$p_elig_for_proj,
                   poly_k = 3)
        
        var_beta_spline_2 <-
          var_beta(df_IF = df_IF,
                   beta = beta_spline_2$x,
                   trial_range = t_range,
                   weights = w,
                   form = 'spline',
                   ate_vec =  df_fold$chi_DR_for_proj,
                   p_elig = df_fold$p_elig_for_proj,
                   spline_knots = sk2)
        
        var_beta_spline_3 <-
          var_beta(df_IF = df_IF,
                   beta = beta_spline_3$x,
                   trial_range = t_range,
                   weights = w,
                   form = 'spline',
                   ate_vec =  df_fold$chi_DR_for_proj,
                   p_elig = df_fold$p_elig_for_proj,
                   spline_knots = sk3)
        
        ### Compute induced PSI
        df_projections <- 
          df_fold %>% 
          mutate('psi_constant' = beta_constant$x,
                 'psi_linear' = as.vector(X_linear %*% beta_linear$x),
                 'psi_cubic' = as.vector(X_cubic %*% beta_cubic$x),
                 'psi_spline_2' = as.vector(X_s2 %*% beta_spline_2$x),
                 'psi_spline_3' = as.vector(X_s3 %*% beta_spline_3$x),
                 'sd_constant' = sqrt(as.vector(var_beta_constant)),
                 'sd_linear' = compute_SD(X_linear, var_beta_linear),
                 'sd_cubic' = compute_SD(X_cubic, var_beta_cubic),
                 'sd_spline_2' = compute_SD(X_s2, var_beta_spline_2),
                 'sd_spline_3' = compute_SD(X_s3, var_beta_spline_3))  %>% 
          pivot_longer(cols = c(contains('psi_'), contains('sd_')),
                       names_to = c('.value', 'projection'),
                       names_pattern = "^([^_]+)_(.+)$") %>%
          rename('psi_hat' = psi) %>%
          mutate('weights' = weights)
        
        proj_summary <- 
          df_projections %>% 
          group_by(projection) %>% 
          summarise('weighted_sd' = weighted.mean(sd, p_elig))
        
        
        ### Evaluate Loss Function
        fold_mse <- 
          df_projections %>% 
          group_by(projection) %>% 
          summarise('mse' = sum( (psi_hat - chi_DR)^2 ),
                    'wELIG_mse' = sum( p_elig * (psi_hat - chi_DR)^2 )) %>% 
          mutate('outcome' = outcome_,
                 'procedures' = procedures_,
                 'weights' = weights)
        
        fold_risk <-
          df_IF %>%
          filter(fold_id == fold) %>%
          select(subject_id, trial_id, chi_IF, eligible) %>%
          inner_join(df_projections, by = 'trial_id', relationship = 'many-to-many') %>%
          group_by(subject_id, projection) %>%
          summarise('risk_IF' = sum(w * (-2 * psi_hat * eligible/p_elig * chi_IF)),
                    'risk_IF_adj' = sum(w * (psi_hat^2 - 2 * psi_hat * eligible/p_elig * chi_IF))) %>%
          group_by(projection) %>%
          summarise('risk' = mean(risk_IF),
                    'risk_adj' = mean(risk_IF_adj),
                    'n_subjects' = n()) %>%
          inner_join(proj_summary, by = 'projection') %>% 
          mutate('outcome' = outcome_,
                 'procedures' = procedures_,
                 'weights' = weights)
        
        df_risk <-
          df_risk %>%
          bind_rows(fold_risk)
        
        df_mse <- 
          df_mse %>% 
          bind_rows(fold_mse)
        
      }
    }
  }
}
write_csv(df_risk, glue('{data_dir}/tv_effects/weight_projection_risk_byfold.csv'))
write_csv(df_mse, glue('{data_dir}/tv_effects/weight_projection_mse_byfold.csv'))

df_risk <- 
  df_risk %>% 
  group_by(outcome, procedures, weights, projection) %>% 
  summarise('risk' = mean(risk),
            'risk_adj' = mean(risk_adj),
            'weighted_sd' = weighted.mean(weighted_sd, sqrt(n_subjects - 1))/sqrt(n()),
            'weighted_sd_adj' = 1/2 * weighted_sd) %>% 
  ungroup() 

best_projections <- 
  df_risk %>% 
  group_by(outcome, procedures, weights) %>% 
  filter(risk == min(risk)) %>% 
  ungroup() 

best_projections_adj <- 
  df_risk %>% 
  group_by(outcome, procedures, weights) %>% 
  filter(risk_adj == min(risk_adj)) %>% 
  ungroup()

df_mse <- 
  df_mse %>% 
  group_by(outcome, procedures, weights, projection) %>% 
  summarise('mse' = mean(mse),
            'wELIG_mse' = mean(wELIG_mse))

best_projections_mse <- 
  df_mse %>% 
  group_by(outcome, procedures, weights) %>% 
  filter(mse == min(mse)) %>% 
  ungroup() 

best_projections_wmse <- 
  df_mse %>% 
  group_by(outcome, procedures, weights) %>% 
  filter(wELIG_mse == min(wELIG_mse)) %>% 
  ungroup()



write_csv(df_risk, glue('{data_dir}/tv_effects/weight_projection_risk.csv'))
write_csv(df_mse, glue('{data_dir}/tv_effects/weight_projection_mse.csv'))
write_csv(best_projections, glue('{data_dir}/tv_effects/best_projections.csv'))
write_csv(best_projections_adj, glue('{data_dir}/tv_effects/best_projections_adj.csv'))
write_csv(best_projections_mse, glue('{data_dir}/tv_effects/best_projections_mse.csv'))
write_csv(best_projections_wmse, glue('{data_dir}/tv_effects/best_projections_wmse.csv'))





### Projection Step (on All Data) + Compute Variance
### Load in Results
df_results <- read_csv(glue('{data_dir}/tv_effects/weight_DR_results.csv'))
df_importance <- read_csv(glue('{data_dir}/tv_effects/weight_DR_importance.csv'))
df_projections <- NULL
for(outcome_ in c('delta_6mo', 'delta_1yr', 'delta_2yr', 'delta_3yr')) {
  for(procedures_ in c('AGB/RYGB/SLEEVE', 'RYGB', 'SLEEVE', 'RYGB/SLEEVE')) {
    cat(outcome_, procedures_, '\n')

    ### Subject Specific IF Contributions
    df_IF <- read_parquet(glue('{data_dir}/tv_effects/IF_contributions/subject_IF_{outcome_}_{gsub("/", "-", procedures_)}.parquet'))

    if(procedures_ == 'AGB/RYGB/SLEEVE') {
      t_range <- 1:84
      sk2 <- c(24, 48)
      sk3 <- c(12, 36, 60)
    } else if(procedures_ == 'RYGB') {
      t_range <- 1:84
      sk2 <- c(24, 48)
      sk3 <- c(12, 36, 60)
    } else if(procedures_ == 'SLEEVE'){
      t_range <- 30:84
      sk2 <- c(48, 72)
      sk3 <- c(40, 56, 72)
    }  else if(procedures_ == 'RYGB/SLEEVE'){
      t_range <- 30:84
      sk2 <- c(48, 72)
      sk3 <- c(40, 56, 72)
    }

    df_info <-
      df_results %>%
      filter(outcome == outcome_) %>%
      filter(procedures == procedures_)

    ### w(m) options
    ### unifom w(m) = 1
    ### inv_elig = 1/P(E[m] = 1) (eg. to cancel out P(E[m] = 1) in estimating eqn
    for(weights in c('uniform')) {
      if(weights == 'uniform') {
        w <- rep(1, length(t_range))
      }

      ### Constant projection
      cat('Constant Projection\n')
      beta_constant <-
        nleqslv(x = c(0),
                fn = estimating_eqn,
                trial_range = t_range,
                weights = w,
                form = 'constant',
                ate_vec = df_info$chi_DR,
                p_elig = df_info$p_elig)

      var_beta_constant <-
        var_beta(df_IF = df_IF,
                 beta = beta_constant$x,
                 trial_range = t_range,
                 weights = w,
                 form = 'constant',
                 ate_vec = df_info$chi_DR,
                 p_elig = df_info$p_elig)

      ### Linear Projection
      cat('Linear Projection\n')
      beta_linear <-
        nleqslv(x = c(0,0),
                fn = estimating_eqn,
                trial_range = t_range,
                weights = w,
                form = 'polynomial',
                ate_vec = df_info$chi_DR,
                p_elig = df_info$p_elig,
                poly_k = 1)

      var_beta_linear <-
        var_beta(df_IF = df_IF,
                 beta = beta_linear$x,
                 trial_range = t_range,
                 weights = w,
                 form = 'polynomial',
                 ate_vec = df_info$chi_DR,
                 p_elig = df_info$p_elig,
                 poly_k = 1)


      ### Cubic Projection
      cat('Cubic Projection\n')
      cubic_start <-
        lm(chi_DR ~ trial_id + I(trial_id^2) + I(trial_id^3),
           data = df_info,
           weight = p_elig)

      beta_cubic <-
        nleqslv(x = coef(cubic_start),
                fn = estimating_eqn,
                trial_range = t_range,
                weights = w,
                form = 'polynomial',
                ate_vec = df_info$chi_DR,
                p_elig = df_info$p_elig,
                poly_k = 3)

      var_beta_cubic <-
        var_beta(df_IF = df_IF,
                 beta = beta_cubic$x,
                 trial_range = t_range,
                 weights = w,
                 form = 'polynomial',
                 ate_vec = df_info$chi_DR,
                 p_elig = df_info$p_elig,
                 poly_k = 3)

      ### Spline Projections
      cat('Spline Projection\n')
      beta_spline_2 <-
        nleqslv(x = c(0,0,0,0),
                fn = estimating_eqn,
                trial_range = t_range,
                weights = w,
                form = 'spline',
                ate_vec = df_info$chi_DR,
                p_elig = df_info$p_elig,
                spline_knots = sk2)

      var_beta_spline_2 <-
        var_beta(df_IF = df_IF,
                 beta = beta_spline_2$x,
                 trial_range = t_range,
                 weights = w,
                 form = 'spline',
                 ate_vec = df_info$chi_DR,
                 p_elig = df_info$p_elig,
                 spline_knots = sk2)

      beta_spline_3 <-
        nleqslv(x = c(0,0,0,0,0),
                fn = estimating_eqn,
                trial_range = t_range,
                weights = w,
                form = 'spline',
                ate_vec = df_info$chi_DR,
                p_elig = df_info$p_elig,
                spline_knots = sk3)

      var_beta_spline_3 <-
        var_beta(df_IF = df_IF,
                 beta = beta_spline_3$x,
                 trial_range = t_range,
                 weights = w,
                 form = 'spline',
                 ate_vec = df_info$chi_DR,
                 p_elig = df_info$p_elig,
                 spline_knots = sk3)

      basis2 <- ns(x = t_range, knots = sk2)
      basis3 <- ns(x = t_range, knots = sk3)

      ### Design Matricies
      X_linear <- cbind(1, t_range)
      X_cubic <- cbind(1, t_range, t_range^2, t_range^3)
      X_s2 <- cbind(1, basis2)
      X_s3 <- cbind(1, basis3)

      ### Compute induced PSI
      df_tmp <-
        df_info %>%
        mutate('psi_constant' = beta_constant$x,
               'psi_linear' = as.vector(X_linear %*% beta_linear$x),
               'psi_cubic' = as.vector(X_cubic %*% beta_cubic$x),
               'psi_spline_2' = as.vector(X_s2 %*% beta_spline_2$x),
               'psi_spline_3' = as.vector(X_s3 %*% beta_spline_3$x),
               'sd_constant' = sqrt(as.vector(var_beta_constant)),
               'sd_linear' = compute_SD(X_linear, var_beta_linear),
               'sd_cubic' = compute_SD(X_cubic, var_beta_cubic),
               'sd_spline_2' = compute_SD(X_s2, var_beta_spline_2),
               'sd_spline_3' = compute_SD(X_s3, var_beta_spline_3)) %>%

        pivot_longer(cols = c(contains('psi_'), contains('sd_')),
                     names_to = c('.value', 'projection'),
                     names_pattern = "^([^_]+)_(.+)$") %>%
        rename('psi_hat' = psi) %>%
        mutate('weights' = weights)

      df_projections <- bind_rows(df_projections, df_tmp)

      cat('\n--------\n')
    }
  }
}

write_csv(df_projections, glue('{data_dir}/tv_effects/weight_projections.csv'))
