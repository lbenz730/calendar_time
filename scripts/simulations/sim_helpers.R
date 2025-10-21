### Conduct Projection Analysis
EIF_projection <- function(df_IF) {
  t_range <- min(df_IF$trial_id):max(df_IF$trial_id)
  
  ### Spline Knots
  sk2 <- c(12, 24)
  sk3 <- c(9, 18, 27)
  
  ### Compute values of CHI for projection and evaluation
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
  
  ### Evaluate Loss Functions
  df_risk <- NULL
  df_sse <- NULL
  
  for(fold in 1:max(df_chi$fold_id)) {
    df_fold <- 
      df_chi %>% 
      filter(fold_id == fold)
    
    for(weights in c('uniform')) {
      if(weights == 'uniform') {
        w <- rep(1, length(t_range))
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
      
      
      
      ### Design Matricies
      basis2 <- ns(x = t_range, knots = sk2)
      basis3 <- ns(x = t_range, knots = sk3)
      X_linear <- cbind(1, t_range)
      X_cubic <- cbind(1, t_range, t_range^2, t_range^3)
      X_s2 <- cbind(1, basis2)
      X_s3 <- cbind(1, basis3)
      
      
      
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
      fold_sse <- 
        df_projections %>% 
        group_by(projection) %>% 
        summarise('mse' = sum( (psi_hat - chi_DR)^2 ),
                  'wELIG_mse' = sum( p_elig * (psi_hat - chi_DR)^2 )) %>% 
        mutate('weights' = weights)
      
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
        mutate('weights' = weights)
      
      df_risk <-
        df_risk %>%
        bind_rows(fold_risk)
      
      df_sse <- 
        df_sse %>% 
        bind_rows(fold_sse)
      
    }
  }
  
  df_risk <- 
    df_risk %>% 
    group_by(weights, projection) %>% 
    summarise('risk' = mean(risk),
              'risk_adj' = mean(risk_adj),
              'weighted_sd' = weighted.mean(weighted_sd, sqrt(n_subjects - 1))/sqrt(n()),
              'weighted_sd_adj' = 1/2 * weighted_sd) %>% 
    ungroup() 
  
  df_sse <- 
    df_sse %>% 
    group_by(weights, projection) %>% 
    summarise('mse' = mean(mse),
              'wELIG_mse' = mean(wELIG_mse)) %>% 
    ungroup()
  
  df_loss <- 
    df_risk %>% 
    inner_join(df_sse, by = c('weights', 'projection'))
  
  
  return(df_loss)
}


### Summary metrics for standardization matrix 
sm_summary_metrics <- function(standardization_matrix, threshold) {
  ### Non-parametric Ratio Statistic
  df_sigma_m <- 
    standardization_matrix %>% 
    pivot_longer(cols = contains('chi'),
                 names_to = 'estimator',
                 values_to = 'chi_hat') %>% 
    group_by(baseline_trial_id, estimator) %>% 
    mutate('mj_diff' = chi_hat - chi_hat[trial_id == baseline_trial_id],
           'mj_diff_threshold' = case_when(abs(mj_diff) >= threshold ~ mj_diff,
                                           T ~ 0)) %>% 
    summarise('sigma2_m' = 1/(n()-1) * sum(mj_diff^2),
              'sigma2_m_threshold' = 1/(n()-1) * sum(mj_diff_threshold^2)) %>% 
    ungroup()
  
  df_gamma_m <- 
    standardization_matrix %>% 
    pivot_longer(cols = contains('chi'),
                 names_to = 'estimator',
                 values_to = 'chi_hat') %>% 
    group_by(trial_id, estimator) %>% 
    mutate('mj_diff' = chi_hat - chi_hat[trial_id == baseline_trial_id],
           'mj_diff_threshold' = case_when(abs(mj_diff) >= threshold ~ mj_diff,
                                           T ~ 0)) %>% 
    summarise('gamma2_m' = 1/(n()-1) * sum(mj_diff^2),
              'gamma2_m_threshold' = 1/(n()-1) * sum(mj_diff_threshold^2)) %>% 
    ungroup()
  
  df_sigma <- 
    df_gamma_m %>% 
    inner_join(df_sigma_m, by = c('trial_id' = 'baseline_trial_id', 'estimator')) %>% 
    mutate('sigma_ratio' = sqrt(sigma2_m)/(sqrt(sigma2_m) + sqrt(gamma2_m)),
           'sigma2_ratio' = sigma2_m/(sigma2_m + gamma2_m),
           'sigma2_ratio_threshold' = sigma2_m_threshold/(sigma2_m_threshold + gamma2_m_threshold)) %>% 
    mutate('sigma_ratio' = ifelse(sigma2_m == 0, 0, sigma_ratio),
           'sigma2_ratio' = ifelse(sigma2_m == 0, 0, sigma2_ratio),
           'sigma2_ratio_threshold' = ifelse(sigma2_m_threshold == 0, 0, sigma2_ratio_threshold))
  
  ### Check how much thresholding mattered
  decisions <- 
    df_sigma %>% 
    group_by(estimator) %>% 
    summarise('avg' = mean(sigma2_ratio),
              'avg_thresh' = mean(sigma2_ratio_threshold)) %>% 
    mutate('decision' = ifelse( abs(avg - avg_thresh) <= 0.1, 'No Threshold', 'Threshold'))
  
  df_sigma <- 
    df_sigma %>% 
    inner_join(decisions, by = 'estimator') %>% 
    mutate('gamma2_m' = ifelse(decision == 'No Threshold', gamma2_m, gamma2_m_threshold),
           'sigma2_m' = ifelse(decision == 'No Threshold', sigma2_m, sigma2_m_threshold)) %>% 
    mutate('sigma_ratio' = sqrt(sigma2_m)/(sqrt(sigma2_m) + sqrt(gamma2_m)),
           'sigma2_ratio' = sigma2_m/(sigma2_m + gamma2_m)) %>% 
    mutate('sigma_ratio' = ifelse(sigma2_m == 0, 0, sigma_ratio),
           'sigma2_ratio' = ifelse(sigma2_m == 0, 0, sigma2_ratio)) %>% 
    select(-contains('threshold'), -contains('avg'))
  
  
  
  
  return(df_sigma) 
  
}


### Custom SL prediction wrapper because by default it doesn't use multithreading for SL.ranger prediction 
SL.predict <- function(model_object, X, binomial = F, return_all = F) {
  preds <- matrix(nrow = nrow(X), ncol = length(model_object$fitLibrary))
  which_keep <- which(model_object$coef != 0) ### ensure we only predict base lnr which have none zero model coeff
  
  #### Extract Component Predictions
  for(i in which_keep) {
    lnr <- model_object$fitLibrary[[i]]$object 
    
    if(grepl('SL.ranger', model_object$libraryNames[i])) {
      if(binomial) {
        preds[,i] <- predict(lnr, data = X)$predictions[,'1']
      } else {
        preds[,i] <- predict(lnr, data = X)$predictions
      }
    } else {
      preds[,i] <- predict(lnr, newdata = X, type = 'response')
    }
  }
  
  if(return_all) {
    return(preds[,which_keep]) 
  }
  
  ### Multiply Component Preds by Model Weights
  if(length(which_keep) > 1) {
    pred_vec <- as.vector(preds[,which_keep] %*% model_object$coef[which_keep])
  } else {
    pred_vec <- as.vector(preds[,which_keep] * model_object$coef[which_keep])
  }
  return(pred_vec)
}

one_hot_encode <- function(df) {
  df <- 
    model.matrix(~., df) %>% 
    as_tibble() %>% 
    select(-1) ### Remove Intercept
  
  return(df)
}

clean_test_df <- function(test, train) {
  ### Select Columns in Train
  test_OHE <- 
    one_hot_encode(bind_rows(test, train)) %>% 
    select(all_of(names(one_hot_encode(train)))) %>% 
    slice(1:nrow(test))
  
  return(test_OHE) 
}


### Fast Function to compute mean sigma/gamma ratio
fast_ratio <- function(vec, M, threshold) {
  standardization_matrix <- matrix(vec, nrow = M, ncol = M, byrow = T)
  
  chi_vec <- diag(standardization_matrix)
  row_null <- outer(chi_vec, rep(1, M))
  col_null <- t(outer(chi_vec, rep(1, M)))
  row_diff <- standardization_matrix - row_null
  col_diff <- standardization_matrix -  col_null
  sigma2 <- rowSums( (row_diff)^2 )/(M-1) 
  gamma2 <- colSums( (col_diff)^2 )/(M-1) 
  
  ### Check for degenerate case against threshold
  row_diff_t <- row_diff
  row_diff_t[abs(row_diff) < threshold] <- 0
  col_diff_t <- col_diff
  col_diff_t[abs(col_diff) < threshold] <- 0
  
  sigma2_t <- rowSums( (row_diff_t)^2 )/(M-1) 
  gamma2_t <- colSums( (col_diff_t)^2 )/(M-1) 
  
  ratio <- mean(sigma2/(sigma2 + gamma2 + 1e-100))
  ratio_t <- mean(sigma2_t/(sigma2_t + gamma2_t + 1e-100))
  
  if(abs(ratio_t - ratio) > 0.1) {
    ratio <- ratio_t
  }
  
  return(ratio)
}



boot_ratio_ci <- function(standardize_IF, standardization_df, threshold, M, n_boot) {
  
  standardization_matrix <- matrix(standardization_df$chi_cross, nrow = M, ncol = M, byrow = T)
  
  chi_vec <- diag(standardization_matrix)
  row_null <- outer(chi_vec, rep(1, M))
  col_null <- t(outer(chi_vec, rep(1, M)))
  row_diff <- standardization_matrix - row_null
  col_diff <- standardization_matrix -  col_null
  sigma2 <- rowSums( (row_diff)^2 )/(M-1) 
  gamma2 <- colSums( (col_diff)^2 )/(M-1) 
  
  ### Check for degenerate case against threshold
  row_diff_t <- row_diff
  row_diff_t[abs(row_diff) < threshold] <- 0
  col_diff_t <- col_diff
  col_diff_t[abs(col_diff) < threshold] <- 0
  
  sigma2_t <- rowSums( (row_diff_t)^2 )/(M-1) 
  gamma2_t <- colSums( (col_diff_t)^2 )/(M-1) 
  
  test_stat <- mean(sigma2/(sigma2 + gamma2 + 1e-100))
  test_stat_t <- mean(sigma2_t/(sigma2_t + gamma2_t + 1e-100))
  
  if(abs(test_stat_t - test_stat) > 0.1) {
    test_stat <- test_stat_t
  }
  
  
  ### Covariance Matrix for Standardization Matrix
  X_IF <- 
    standardize_IF %>% 
    select(subject_id, baseline_trial_id, trial_id, chi_IFj) %>% 
    mutate('baseline_trial_id' = paste0('j', baseline_trial_id),
           'trial_id' = paste0('m', trial_id)) %>% 
    pivot_wider(names_from = c('trial_id', 'baseline_trial_id'),
                values_from = 'chi_IFj') %>% 
    select(-subject_id) %>% 
    as.matrix()
  
  X_IF <- scale(X_IF, center = T, scale = F)   # center columns
  Sigma_IF <- ( crossprod(X_IF) / (nrow(X_IF) - 1) ) * 1/nrow(X_IF)
  
  
  ### Parametric Bootstrap
  resamples <- MASS::mvrnorm(n = n_boot, mu = as.vector(t(standardization_matrix)), Sigma = Sigma_IF) 
  sim_stats <- map_dbl(1:n_boot, ~fast_ratio(vec = resamples[.x,], M, threshold))
  q <- quantile(sim_stats, c(0.025, 0.975))
  
  df_test <- 
    tibble('test_stat' = test_stat,
           'sd_stat' = sd(sim_stats),
           'qlow' = q[1],
           'qhigh' = q[2])
  
  return(df_test)
}
