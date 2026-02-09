
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

## MSE Metrics
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

## Diagnostics:
## Selection proxy: indices of top-k coefficients by absolute magnitude
##
## Note: 
## Ridge does not perform variable selection. 
## This top-k rule is used only as a *stability proxy* to compare which coefficients are emphasized across imputations.
##
## Settings:
## beta_hat  : estimated coefficient vector (length p)
## s_select  : k (number of top coefficients to return)
##
## Returns: sorted integer vector of selected indices

selected_topk <- function(beta_hat, s_select = 10) {
  o <- order(abs(beta_hat), decreasing = TRUE)
  sort(o[seq_len(min(s_select, length(o)))])
}

## Jaccard similarity between two sets of indices
##
## Settings:
## A, B : integer vectors (selected sets)
##
## Returns:
## J(A,B) = |A ∩ B| / |A ∪ B|, NA if union is empty

jaccard <- function(A, B) {
  A <- unique(A); B <- unique(B)
  inter <- length(intersect(A, B))
  uni <- length(union(A, B))
  if (uni == 0) return(NA_real_)
  inter / uni
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
