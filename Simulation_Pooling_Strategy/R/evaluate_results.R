# Simulation_Pooling_Strategy/R/evaluate_predictions.R
# Purpose: Use pooled coefficient vectors to make predictions on the complete test set and compute the corresponding MSE values.

# Predict from a pooled coefficient vector
# beta_vec  : pooled coefficient vector
# test_data : complete test data frame
# Returns a numeric vector of predictions on the test set

predict_from_beta <- function(beta_vec, test_data) {
  
  # Extract predictor matrix from test data
  X_test <- as.matrix(test_data[, colnames(test_data) != "y", drop = FALSE])
  
  # Extract intercept and slope coefficients
  intercept <- beta_vec["(Intercept)"]
  beta_x <- beta_vec[setdiff(names(beta_vec), "(Intercept)")]
  
  # Compute predictions
  y_pred <- as.numeric(intercept + X_test %*% beta_x)
  
  return(y_pred)
}

# Compute mean squared error
# y_true : true response vector
# y_pred : predicted response vector
# Returns MSE

compute_mse <- function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

# Evaluate pooled beta vectors on the complete test set
# pooled_beta : list containing pooled coefficient vectors from different methods
# test_data   : complete test data frame
# Returns: A numeric vector of test MSE values for each pooling method

evaluate_pooled_betas <- function(pooled_beta, test_data) {
  y_test <- test_data$y
  mse_out <- sapply(pooled_beta, function(beta_vec) {
    y_pred <- predict_from_beta(beta_vec, test_data)
    compute_mse(y_test, y_pred)
  })
  return(mse_out)
}

# Compute relative MSE ratios using the complete-data benchmark
# mse_pooled   : numeric vector of pooled-method MSE values
# mse_complete : MSE from the complete data
# Returns a numeric vector of relative MSE ratios
compute_relative_mse <- function(mse_pooled, mse_complete) {
  mse_pooled / mse_complete
}