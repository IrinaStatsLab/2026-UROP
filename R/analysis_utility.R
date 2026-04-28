
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

# Compute MSE
mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

# Run MICE with predictive mean matching
run_mice_pmm <- function(df, m, seed) {
  tryCatch({
    ini <- mice::mice(df, m = 1, maxit = 0, printFlag = FALSE)
    method <- ini$method
    pred   <- ini$predictorMatrix
    
    method["y"] <- ""
    method[setdiff(names(df), "y")] <- "pmm"
    pred[, "y"] <- 1
    pred["y", ] <- 0
    
    mice::mice(
      df,
      m = m,
      method = method,
      predictorMatrix = pred,
      printFlag = FALSE,
      seed = seed
    )
  }, error = function(e) NULL)
}

# Select ridge penalty using cross-validation
select_ridge_lambda_cv <- function(X, y, seed = NULL, nfolds = 5) {
  if (!is.null(seed)) set.seed(seed)
  
  cv_fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = as.numeric(y),
    alpha = 0,
    nfolds = nfolds,
    standardize = TRUE
  )
  
  as.numeric(cv_fit$lambda.min)
}

