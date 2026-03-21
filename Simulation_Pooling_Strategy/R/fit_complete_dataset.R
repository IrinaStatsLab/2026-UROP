# Simulation_Pooling_Strategy/R/fit_complete_dataset.R
# Purpose: Fit ridge regression on the complete training data, select the parameter lambda using the complete validation set, 
# and evaluate the MSE on the complete test set.

# Fit ridge regression on complete data and use the validation set to select lambda.
# train_data      : complete training data frame
# validation_data : complete validation data frame
# test_data       : complete test data frame
# lambda_grid     : candidate lambda values
# Returns a list containing:
# lambda              : selected lambda
# fit_path            : fitted ridge object from glmnet
# beta                : coefficient vector at selected lambda
# validation_mse_path : validation MSE for each lambda
# mse_test            : test MSE at selected lambda


fit_complete_benchmark <- function(train_data, validation_data, test_data, lambda_grid) {

  # Extract predictor matrices and response vectors
  X_train <- as.matrix(train_data[, colnames(train_data) != "y", drop = FALSE])
  y_train <- train_data$y
  
  X_val <- as.matrix(validation_data[, colnames(validation_data) != "y", drop = FALSE])
  y_val <- validation_data$y
  
  X_test <- as.matrix(test_data[, colnames(test_data) != "y", drop = FALSE])
  y_test <- test_data$y
  
  # Fit ridge on complete training data
  fit_path <- glmnet::glmnet(x = X_train, y = y_train, alpha = 0, lambda = lambda_grid, standardize = TRUE)
  
  # Predict on validation set for all candidate lambdas
  pred_val <- predict(fit_path, newx = X_val)
  
  # Compute validation MSE path
  validation_mse_path <- apply(pred_val, 2, function(pred) {
    mean((y_val - pred)^2)
  })
  
  # Choose lambda minimizing validation MSE
  best_idx <- which.min(validation_mse_path)
  lambda_best <- lambda_grid[best_idx]
  
  # Extract coefficients at selected lambda
  beta_best <- as.matrix(coef(fit_path, s = lambda_best))
  
  # Predict on test set using selected lambda
  pred_test <- as.numeric(predict(fit_path, newx = X_test, s = lambda_best))
  mse_test <- mean((y_test - pred_test)^2)
  
  return(list(
    lambda = lambda_best,
    fit_path = fit_path,
    beta = beta_best,
    validation_mse_path = validation_mse_path,
    mse_test = mse_test
  ))
}