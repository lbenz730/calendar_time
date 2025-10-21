### Projection Function (stores gradient vector as well)
### m = trial index
### trial_range = range of trial indices for analysis
### beta = vector of MSM coefficients
### form = descriptor of the functional form of the MSM
### optional arguments for MSM structure
###     poly_k = degree of polynomial (for polynomial form)
###     spline_knots = knots for natural cubic splines (for spline form)
psi_projection <- function(m, trial_range, beta, form, poly_k = NULL, spline_knots = NULL) {
  ### Constant
  if(form == 'constant') {
    psi <- beta  
    gradient <- 1
    hessian <- 0
  } else if(form == 'polynomial') {
    poly_basis <- m^(0:poly_k) 
    psi <- as.vector(poly_basis %*% beta)
    gradient <- poly_basis
    hessian <- matrix(data = 0, nrow = length(beta), ncol = length(beta))
  } else if(form == 'spline') {
    spline_basis <- c(1, as.vector(ns(x = m, knots = spline_knots, Boundary.knots = range(trial_range))))
    psi <- as.vector(spline_basis %*% beta)
    gradient <- spline_basis
    hessian <- matrix(data = 0, nrow = length(beta), ncol = length(beta))
  }
  return(list('psi' = psi, 'gradient' = gradient, 'hessian' = hessian))
}



### Computes the given value of the summand for the estimating equation for beta vectors
### beta = current value of beta
### M = max trial_index
### weights = weight vector (w) for estimating equation
### form = descriptor of the functional form of the MSM
### ate_vec = vector of estimates of chi(m) 
### p_elig = eligibility fraction
### optional arguments for MSM structure
###     poly_k = degree of polynomial (for polynomial form)
###     spline_knots = knots for natural cubic splines (for spline form)
estimating_eqn <- function(beta, trial_range, weights, form, ate_vec, p_elig, poly_k = NULL, spline_knots = NULL) {
  proj_info <- map(trial_range, ~psi_projection(m = .x, trial_range, beta, form, poly_k, spline_knots))
  psi_vals <- map_dbl(proj_info, pluck, 'psi') 
  gradient_matrix <- do.call(rbind, map(proj_info, pluck, 'gradient'))
  
  ### Summand for estimating equation
  ee_summand <- (weights * p_elig * (ate_vec - psi_vals)) %*% gradient_matrix
  return(ee_summand)
}

### Computes the normalization matrix V necessary to characterize the variance of beta
normalization_V <- function(beta, trial_range, weights, form, ate_vec, p_elig, poly_k = NULL, spline_knots = NULL) {
  proj_info <- map(trial_range, ~psi_projection(m = .x, trial_range, beta, form, poly_k, spline_knots))
  psi_vals <- map_dbl(proj_info, pluck, 'psi') 
  gradient_matrices <- map(proj_info, ~{.x$gradient %*% t(.x$gradient) })
  
  summand_matrices <- 
    map(1:length(trial_range), function(m) {
      weights[m] * p_elig[m] * ( proj_info[[m]]$hessian  * (ate_vec[m] - psi_vals[m])  - gradient_matrices[[m]])
    })
  
  V <- Reduce('+', summand_matrices)
  
  return(V)
}

### Computes "meat" portion of variance estimator for variance of beta
sigma_compute <- function(df_IF, beta, trial_range, weights, form, ate_vec, poly_k = NULL, spline_knots = NULL) {
  proj_info <- map(trial_range, ~psi_projection(m = .x, trial_range, beta, form, poly_k, spline_knots))
  psi_vals <- map_dbl(proj_info, pluck, 'psi') 
  gradient_matrix <- do.call(rbind, map(proj_info, pluck, 'gradient'))
  
  ### Join in Info for efficient Computation
  df_info <- 
    tibble('w' = weights, 
           'm' = trial_range,
           'psi' = psi_vals)
  
  x <- 
    df_IF %>% 
    inner_join(df_info, by = c('trial_id' = 'm')) %>% 
    mutate('contribution' = w * eligible * (chi_IF - psi)) %>% 
    select(subject_id, trial_id, contribution) %>% 
    pivot_wider(names_from = 'trial_id',
                values_from = 'contribution') %>% 
    select(subject_id, order(as.numeric(names(.)[-1])) + 1)
  
  ### Matrix of beta IF contributions
  beta_IF <- as.matrix(x[,-1]) %*% gradient_matrix
  Sigma <- cov(beta_IF)
  
  return(Sigma)
}

var_beta <- function(df_IF, beta, trial_range, weights, form, ate_vec, p_elig, poly_k = NULL, spline_knots = NULL) {
  V <- normalization_V(beta, trial_range, weights, form, ate_vec, p_elig, poly_k, spline_knots)
  Sigma <- sigma_compute(df_IF, beta, trial_range, weights, form, ate_vec, poly_k, spline_knots)
  sigma_beta <- 1/n_distinct(df_IF$subject_id) * solve(V) %*% Sigma %*% solve(V)
  return(sigma_beta)
}


compute_SD <- function(X, Sigma) {
  std_devs <- map_dbl(1:nrow(X), ~sqrt( t(X[.x,]) %*% Sigma %*% X[.x, ]) )
  return(std_devs)
}
