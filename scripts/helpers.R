library(tidyverse) 

### Inverse Logit Function
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

### Logit Function
logit <- function(x) {
  return(log(x/(1-x)))
}

### Quick way to compute logistic regression coefficients when fitting saturated model with one coefficient
glm_quick <- function(Y, A, W) {
  beta0 <- logit(weighted.mean(Y[A == 0], w = W[A == 0])) 
  beta1 <- logit(weighted.mean(Y[A == 1], w = W[A == 1])) - beta0
  return(list('beta0' = beta0, 'beta1' = beta1))
}

### Function to compute the fitted value of a model given a list of coefficients
compute_model <- function(df, beta) {
  
  ### Initialize fitted values to be 0
  preds <- rep(0, nrow(df))
  
  ### Iterate over coefficients
  for(i in 1:length(beta)) {
    coeff <- names(beta)[i]
    
    if(coeff == '(Intercept)') { ### Intercept
      preds <- preds + beta[[coeff]]
    } else if(grepl('\\[', coeff) & grepl(':', coeff)) { ### Level of categorical variable w/ interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeffs[1]))
      variable <- gsub('\\[.*$', '', coeffs[1])
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + df[[ coeffs[2] ]][ix] * beta[[ coeff ]]
    } else if(grepl('\\[', coeff)) { ### Level of categorical variable
      level <- gsub('^.*\\[', '',  gsub('\\]', '', coeff))
      variable <- gsub('\\[.*$', '', coeff)
      ix <- df[[variable]] == level
      preds[ix] <- preds[ix] + beta[[coeff]]
    } else if(grepl(':', coeff)) { ### Interaction
      coeffs <- unlist(strsplit(coeff, ':'))
      preds <- preds + df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * beta[[coeff]]
    } else if(grepl('\\^2\\)$', coeff)) { ### quadratic 
      variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
      preds <- preds + df[[ variable ]]^2 * beta[[coeff]]
    } else if(grepl('^I\\(exp\\(', coeff)) { ### exponential
      variable <- gsub('\\)\\)$', '', gsub('^I\\(exp\\(', '', coeff))
      preds <- preds + exp(df[[ variable ]]) * beta[[coeff]]
    } else { ### Linear Term 
      preds <- preds + df[[coeff]] * beta[[coeff]]
    }
  }
  
  return(preds)
}

### Function to compute the matrix of fitted group probabilities for a 
### multinomial given a list of coefficients
compute_multinomial_model <- function(df, beta) {
  ### Initialize to logit of cumumlative probability of being in each group
  probs <- matrix(beta$intercepts, nrow = nrow(df), ncol = length(beta$intercepts), byrow = T)
  
  ### Add in rest of covariate effects
  for(i in 1:length(beta)) {
    coeff <- names(beta)[i]
    
    ### Note that positive beta --> subtraction 
    ### See: https://www.bookdown.org/rwnahhas/RMPH/blr-ordinal.html
    if(!(coeff %in% c('intercepts', '(Intercept)'))) {
      if(grepl('\\[', coeff) & grepl(':', coeff)) { ### Level of categorical variable w/ interaction
        coeffs <- unlist(strsplit(coeff, ':'))
        level <- gsub('^.*\\[', '',  gsub('\\]', '', coeffs[1]))
        variable <- gsub('\\[.*$', '', coeffs[1])
        ix <- df[[variable]] == level
        probs[ix,] <- probs[ix] + df[[ coeffs[2] ]][ix] * beta[[ coeff ]]
      } else if(grepl(':', coeff)) { ### Interaction
        coeffs <- unlist(strsplit(coeff, ':'))
        probs <- probs - df[[ coeffs[1] ]] * df[[ coeffs[2] ]] * beta[[coeff]]
      }  else if(grepl('\\[', coeff)) { ### Level of categorical variable
        level <- gsub('^.*\\[', '',  gsub('\\]', '', coeff))
        variable <- gsub('\\[.*$', '', coeff)
        ix <- df[[variable]] == level
        probs[ix,] <- probs[ix,] - beta[[coeff]]
      } else if(grepl('\\^2\\)$', coeff)) { ### quadratic 
        variable <- gsub('\\^2\\)$', '', gsub('^I\\(', '', coeff))
        probs <- probs - df[[ variable ]]^2 * beta[[coeff]]
      } else if(grepl('^I\\(exp\\(', coeff)) { ### exponential
        variable <- gsub('\\)\\)$', '', gsub('^I\\(exp\\(', '', coeff))
        probs <- probs - exp(df[[ variable ]]) * beta[[coeff]]
      } else { ### Linear Term 
        probs <- probs - df[[coeff]] * beta[[coeff]]
      }
    }
  }
  
  return(probs)
}

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 12),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))


### Clean names so we know what variables to keep in our simulations
clean_names <- function(params, df_names) {
  x <- unique(gsub('^.*\\.', '', names(unlist(params))))
  x <- gsub('\\[.*', '', unlist(strsplit(x, ':')))
  x <- c(x, unique(unlist(map(params$truth_model, as.character))))
  x <- x[x %in% df_names]
  return(x)
}

### Function to replace null values w/ useful string
replace_null <- function(x, replace_string = NA) {
  x[is.null(x)] <- replace_string
  return(x)
}


winsorize <- function(x, q) {
  ### Get upper and lower quantiles of the data
  q_data <- quantile(x[!is.na(x)], q)
  lower <- q_data[1]
  upper <- q_data[2]
  
  ### Replace extreme values by quantiles
  x[which(x < lower)] <- lower
  x[which(x > upper)] <- upper
  
  return(x)
  
}


clean_test_matrix <- function(test, train) {
  ### Clean Names
  colnames(test) <-  gsub(':', '_', colnames(test))
  
  ### Remove columns in test that aren't in train
  test <- test[,intersect(colnames(test), colnames(train))]
  
  ### Add columns in test that are in train but not test
  diff_names <- setdiff(colnames(train), colnames(test))
  if(length(diff_names) > 0) {
    M <- 
      matrix(data = 0, 
             nrow = nrow(test),
             ncol = length(diff_names))
    colnames(M) <- diff_names
    test <- cbind(test, M)
    test <- test[,colnames(train)]
    
  }
  return(test)
}

### Convert model to structure of list 
model_list <- function(model) {
  broom::tidy(model) %>% 
    mutate('output' = paste0("'", term, "' = ", sprintf('%0.3f', estimate))) %>% 
    pull(output) %>% 
    cat(., sep = ',\n')
}


ar1_Sigma <- function(n, rho, sigma) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  return(sigma^2 * rho^exponent)
}


### CKD EPI Equation for converting serum creatinine into estimated glucose filtration rate
ckd_epi <- function(scr, age, sex) {
  scr <- pmax(0.1, scr)
  alpha <- ifelse(sex == 'M', -0.302, -0.241)
  kappa <- ifelse(sex == 'M', 0.9, 0.7)
  
  eGFR <- 142 * pmin(scr/kappa, 1)^alpha * pmax(scr/kappa, 1)^(-1.2) *  0.9938^age * ifelse(sex == 'M', 1, 1.012)
  return(eGFR)
}

### Emulate ggplot defaults https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### Work Around for gtsave on the cluster
cluster_gtsave <- function(table, file_path, keep_html = F, ...) {
  file <- gsub('\\.html', '', gsub('\\.png', '', file_path))
  gtsave(table, paste0(file, '.html'))
  webshot::webshot(paste0(file, '.html'), paste0(file, '.png'), ...)
  if(!keep_html) {
    file.remove(paste0(file, '.html'))
  }
}
