
# simulation/R/analysis_utility.R
# Contains: ridge regression fit, MSE metrics, selection stability diagnostics, Mice wrapper

## Fit Ridge Regression Model with a Fixed Penalty (lambda)
## Settings:
## X_fit        : n x p numeric matrix 
## y_fit        : response vector (length n)
## lambda_fixed : fixed ridge penalty (scalar > 0)
##
## Returns: glmnet model object (with alpha = 0 for ridge)

fit_ridge_fixed_lambda <- function(X_fit, y_fit, lambda_fixed) {
  glmnet::glmnet(
    x = X_fit, y = y_fit, alpha = 0,
    lambda = lambda_fixed, standardize = TRUE
  )
}

## Predict from a fitted glmnet ridge model at its fitted lambda
##
## Settings:
## fit    : glmnet model object returned by fit_ridge_fixed_lambda()
## X_eval : n_eval x p matrix for evaluation
##
## Returns: numeric vector of predictions (length n_eval)

predict_ridge <- function(fit, X_eval) {
  as.numeric(predict(fit, newx = X_eval, s = fit$lambda))
}

## Extract coefficients from ridge fit
coef_ridge <- function(fit) {
  beta <- as.matrix(stats::coef(fit, s = fit$lambda))
  beta <- beta[setdiff(rownames(beta), "(Intercept)"), , drop = FALSE]
  out <- as.numeric(beta[, 1])
  names(out) <- rownames(beta)
  out
}

## MSE Metrics
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

## Diagnostics:
## Permutation-threshold feature selection for RIDGE
#  Parameters: 
## X: n x p numeric matrix or data.frame of predictors
## y: numeric response vector of length n
## lambda_fixed: fixed ridge penalty parameter λ
## perm_m: integer > 0; number of permutations used to estimate the null distribution of max |beta|
## sel_alpha: numeric in (0,1); quantile level for threshold selection

ridge_select_perm_threshold <- function(X, y, lambda_fixed,
                                        perm_m = 30, sel_alpha = 0.95,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  lambda_fixed <- as.numeric(lambda_fixed)[1]
  
  # Fit on observed y
  fit_obs  <- fit_ridge_fixed_lambda(X_fit = X, y_fit = y, lambda_fixed = lambda_fixed)
  beta_obs <- coef_ridge(fit_obs)
  
  # Null distribution: permute y, refit ridge, record max |beta|
  null_max <- numeric(perm_m)
  for (m in seq_len(perm_m)) {
    y_perm <- sample(y, replace = FALSE)
    fit_null  <- fit_ridge_fixed_lambda(X_fit = X, y_fit = y_perm, lambda_fixed = lambda_fixed)
    beta_null <- coef_ridge(fit_null)
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
