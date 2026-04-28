
# R/analysis_utility.R
# Utility functions for ridge regression, prediction, and MSE evaluation

# Fit ridge regression with a fixed penalty (Lambda)
fit_ridge_fixed_lambda <- function(X_fit, y_fit, lambda_fixed) {
  glmnet::glmnet(
    x = X_fit, y = y_fit, alpha = 0,
    lambda = lambda_fixed, standardize = TRUE)
}

# Predict responses from a fitted ridge model
predict_ridge <- function(fit, X_eval) {
  as.numeric(predict(fit, newx = X_eval, s = fit$lambda))
}

# Extract coefficient vector (excluding intercept)
coef_ridge <- function(fit) {
  beta <- as.matrix(stats::coef(fit, s = fit$lambda))
  beta <- beta[setdiff(rownames(beta), "(Intercept)"), , drop = FALSE]
  out <- as.numeric(beta[, 1])
  names(out) <- rownames(beta)
  out
}

# Extract intercept from ridge model
intercept_ridge <- function(fit) {
  beta_all <- as.matrix(stats::coef(fit, s = fit$lambda))
  as.numeric(beta_all["(Intercept)", 1])
}

# Predict using pooled coefficients (used after pooling across imputations)
predict_from_pooled_beta <- function(beta, intercept, X_eval) {
  as.numeric(intercept + X_eval %*% beta)
}