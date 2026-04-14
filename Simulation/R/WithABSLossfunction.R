
# simulation/R/analysis_utility.R
# Contains: lasso regression fit, MSE metrics, selection stability diagnostics, Mice wrapper

## Fit lasso Regression Model with a Fixed Penalty (lambda)
## Settings:
## X_fit        : n x p numeric matrix 
## y_fit        : response vector (length n)
## lambda_fixed : fixed lasso penalty (scalar > 0)
##
## Returns: glmnet model object (with alpha = 1 for Lasso)

fit_lasso_fixed_lambda <- function(X_fit, y_fit, lambda_fixed) {
  glmnet::glmnet(
    x = X_fit, y = y_fit, alpha = 1,
    lambda = lambda_fixed, standardize = TRUE)
}

## Predict from a fitted glmnet lasso model at its fitted lambda
##
## Settings:
## fit    : glmnet model object returned by fit_lasso_fixed_lambda()
## X_eval : n_eval x p matrix for evaluation
##
## Returns: numeric vector of predictions (length n_eval)

predict_lasso <- function(fit, X_eval) {
  as.numeric(predict(fit, newx = X_eval, s = fit$lambda))
}

## Extract coefficients from lasso fit
coef_lasso <- function(fit) {
  beta <- as.matrix(stats::coef(fit, s = fit$lambda))
  beta <- beta[setdiff(rownames(beta), "(Intercept)"), , drop = FALSE]
  out <- as.numeric(beta[, 1])
  names(out) <- rownames(beta)
  out
}

##Var and SD estimation by regular bootstrapping
bootstrap_lasso_mse <- function(X, y, B,lambda){
  
  n <- nrow(X)
  mse_vec <- numeric(B)
  
  for(b in 1:B){
    
    # bootstrap sample
    idx <- sample(1:n, n, replace = TRUE)
    
    X_boot <- X[idx, ]
    y_boot <- y[idx]
    
    # out-of-bag samples
    oob <- setdiff(1:n, unique(idx))
    
    if(length(oob) == 0) next
    
    X_test <- X[oob, ]
    y_test <- y[oob]
    
    # fit lasso
    fit <- fit_lasso_fixed_lambda(X_fit = X_boot, y_fit = y_boot,lambda_fixed = lambda )
    
    # prediction
    pred <- predict_lasso(fit = fit, X_eval = X_test)
    
    # compute MSE
    mse_vec[b] <- mean((y_test - pred)^2)
    
  }
  
  mse_vec <- mse_vec[mse_vec != 0]
  
  list(
    mean_MSE = mean(mse_vec),
    var_MSE = var(mse_vec),
    sd_MSE = sd(mse_vec),
    mse_samples = mse_vec
  )
}






## MSE Metrics
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

## Diagnostics:
## Permutation-threshold feature selection for lasso
#  Parameters: 
## X: n x p numeric matrix or data.frame of predictors
## y: numeric response vector of length n
## lambda_fixed: fixed lasso penalty parameter λ
## perm_m: integer > 0; number of permutations used to estimate the null distribution of max |beta|
## sel_alpha: numeric in (0,1); quantile level for threshold selection

lasso_select_perm_threshold <- function(X, y, lambda_fixed,
                                        perm_m = 30, sel_alpha = 0.95,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  lambda_fixed <- as.numeric(lambda_fixed)[1]
  
  # Fit on observed y
  fit_obs  <- fit_lasso_fixed_lambda(X_fit = X, y_fit = y, lambda_fixed = lambda_fixed)
  beta_obs <- coef_lasso(fit_obs)
  
  # Null distribution: permute y, refit lasso, record max |beta|
  null_max <- numeric(perm_m)
  for (m in seq_len(perm_m)) {
    y_perm <- sample(y, replace = FALSE)
    fit_null  <- fit_lasso_fixed_lambda(X_fit = X, y_fit = y_perm, lambda_fixed = lambda_fixed)
    beta_null <- coef_lasso(fit_null)
    null_max[m] <- max(abs(beta_null), na.rm = TRUE)
  }
  
  threshold <- as.numeric(stats::quantile(null_max, probs = sel_alpha, names = FALSE, type = 7))
  selected  <- abs(beta_obs) > threshold
  
  list(selected = selected, threshold = threshold, beta = beta_obs, null_max = null_max)
}

## Wrapper for MICE using PMM(Predictive Mean Matching)
##
## Purpose: Run multiple imputation while preventing the entire simulation from crashing when mice() fails 
## (e.g., due to extreme missingness).
##
## Settings:
## df   : data.frame containing predictors (with NAs) and any fully observed vars (e.g., y)
## m    : number of imputations
## seed : Seed for reproducibility
##
## Returns: mice::mids object if successful; otherwise NULL

run_mice_pmm <- function(df, m, seed) {
  tryCatch(
    mice::mice(df, m = m, method = "pmm", printFlag = FALSE, seed = seed),
    error = function(e) NULL
  )
}

library(dplyr)
library(ggplot2)
library(readr)

df1 <- MAR_and_DENSE %>%
  mutate(mechanism = "MAR", beta_type = "Dense")

df2 <- MAR_and_Sparse %>%
  mutate(mechanism = "MAR", beta_type = "Sparse")

df3 <- MCAR_and_DENSE %>%
  mutate(mechanism = "MCAR", beta_type = "Dense")

df4 <- MCAR_and_Sparse %>%
  mutate(mechanism = "MCAR", beta_type = "Sparse")

df <- bind_rows(df1, df2, df3, df4)

ggplot(df, aes(x = factor(pi), y = ratio_R_mean)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_grid(beta_type ~ mechanism) +
  labs(
    x = "Missingness proportion (π)",
    y = "R = MSE_pooled / MSE_full"
  ) +
  theme_minimal(base_size = 14)

  