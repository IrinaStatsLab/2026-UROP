# Simulation_Pooling_Strategy/R/fit_imputed_datasets.R
# Purpose:
# For each imputed training dataset:
# 1. fit a ridge regression over the candidate lambda grid
# 2. use the complete validation set to select lambda_j
# 3. extract the coefficient vector beta_j at the selected lambda_j

# Fit ridge regression on each imputed training dataset and use the same complete validation set to select lambda.
# imputed_list    : list of completed training datasets
# validation_data : complete validation data frame
# lambda_grid     : candidate lambda values

# Returns a list containing:
# lambda_vec           : selected lambda for each imputed dataset
# beta_mat             : coefficient matrix (rows = coefficients, cols = imputations)
# validation_mse_list  : list of validation MSE 
# fit_list             : list of fitted ridge objects

fit_imputed_datasets <- function(imputed_list, validation_data, lambda_grid) {
  
  # Extract validation predictors and response variable
  X_val <- as.matrix(validation_data[, colnames(validation_data) != "y", drop = FALSE])
  y_val <- validation_data$y
  
  # Initialize the variable list
  n_imp <- length(imputed_list)
  fit_list <- vector("list", n_imp)
  validation_mse_list <- vector("list", n_imp)
  lambda_vec <- numeric(n_imp)
  beta_list <- vector("list", n_imp)
  
  # Loop through each imputed training dataset
  for (j in seq_len(n_imp)) {
    dat_j <- imputed_list[[j]]
    
    # Extract predictors and response for imputed dataset j
    X_train_j <- as.matrix(dat_j[, colnames(dat_j) != "y", drop = FALSE])
    y_train_j <- dat_j$y
    
    # Fit ridge regression path on imputed dataset j
    # alpha = 0 ensures ridge regression
    fit_j <- glmnet::glmnet(x = X_train_j, y = y_train_j, alpha = 0, lambda = lambda_grid, standardize = TRUE)
    
    # Predict on the complete validation set for all candidate lambdas
    pred_val_j <- predict(fit_j, newx = X_val)
    
    # Compute validation MSE path for imputation j
    mse_path_j <- apply(pred_val_j, 2, function(pred) {
      mean((y_val - pred)^2)
    })
    
    # Select lambda_j minimizing validation MSE
    best_idx_j <- which.min(mse_path_j)
    lambda_j <- lambda_grid[best_idx_j]
    
    # Extract coefficient vector at the selected lambda_j
    beta_j <- as.matrix(coef(fit_j, s = lambda_j))
    
    # Store outputs
    fit_list[[j]] <- fit_j
    validation_mse_list[[j]] <- mse_path_j
    lambda_vec[j] <- lambda_j
    beta_list[[j]] <- beta_j
  }
  
  # Combine all coefficient vectors into one matrix
  beta_mat <- do.call(cbind, beta_list)
  colnames(beta_mat) <- paste0("imp_", seq_len(n_imp))
  
  return(list(
    lambda_vec = lambda_vec,
    beta_mat = beta_mat,
    validation_mse_list = validation_mse_list,
    fit_list = fit_list
  ))
}