# Simulation_Pooling_Strategy/R/one_simulation_replicate.R
# Purpose: Run one complete simulation replicate under one scenario.

run_one_simulation <- function(setting) {
  
  # Generate complete data
  complete_obj <- generate_complete_data(n = setting$n,b = setting$b,K = setting$K,
                                         rho = setting$rho, beta_type = setting$beta_type,
                                         s = setting$s, signal = setting$signal, sigma_eps = setting$sigma_eps,
                                         seed = setting$seed_data)
  
  # Split into training / validation / test sets
  split_obj <- split_data(complete_obj = complete_obj, train_prop = setting$train_prop,val_prop = setting$val_prop,seed = setting$seed_split)
  train_data_complete <- split_obj$train
  validation_data_complete <- split_obj$validation
  test_data_complete <- split_obj$test
  
  # Fit ridge benchmark on complete training data
  bench_out <- fit_complete_benchmark(train_data = train_data_complete, validation_data = validation_data_complete, 
                                      test_data = test_data_complete,lambda_grid = setting$lambda_grid)
  mse_complete <- bench_out$mse_test
  
  # Impose missingness on training set only
  train_data_missing <- impose_missingness(data = train_data_complete, mechanism = setting$mechanism,
                                           pi = setting$pi, M_risk = setting$M_risk, a1 = setting$a1, a2 = setting$a2, a3 = setting$a3,
                                           seed = setting$seed_missing)
  
  # Run MICE on incomplete training data
  imp_out <- run_mice_imputation(data_mis = train_data_missing, m = setting$m, maxit = setting$maxit, 
                                 method = setting$imputation_method, seed = setting$seed_mice)
  imputed_list <- imp_out$imputed_list
  
  # Fit ridge on each imputed training dataset
  imputed_out <- fit_imputed_datasets(imputed_list = imputed_list, validation_data = validation_data_complete, lambda_grid = setting$lambda_grid)
  beta_mat <- imputed_out$beta_mat
  
  # Pool coefficient vectors
  pooled_beta <- pool_beta_all(beta_mat)
  
  # Evaluate pooled coefficients on complete test set
  mse_pooled <- evaluate_pooled_betas(pooled_beta = pooled_beta, test_data = test_data_complete)
  
  # Compute relative MSE ratios
  R <- compute_relative_mse(mse_pooled = mse_pooled, mse_complete = mse_complete)
  
  # Compute Coefficient instability across imputations
  # Remove intercept row before computing across-imputation variance
  beta_within_imp_var <- mean(apply(beta_mat[-1, , drop = FALSE], 1, var, na.rm = TRUE))
  
  # Compute the Distance between pooled beta and true beta
  true_beta_full <- c("(Intercept)" = 0, setNames(complete_obj$beta, paste0("X", seq_along(complete_obj$beta))))
  l2_delta_beta_mean <- sqrt(sum((pooled_beta$mean[names(true_beta_full)] - true_beta_full)^2))
  l2_delta_beta_geo <- sqrt(sum((pooled_beta$geo[names(true_beta_full)] - true_beta_full)^2))
  l2_delta_beta_svd <- sqrt(sum((pooled_beta$svd[names(true_beta_full)] - true_beta_full)^2))
  
  # Organize outputs
  results <- data.frame(n = setting$n,b = setting$b,K = setting$K,p = setting$b * setting$K,
                        rho = setting$rho,beta_type = setting$beta_type,s = setting$s,signal = setting$signal,
                        sigma_eps = setting$sigma_eps,mechanism = setting$mechanism,pi = setting$pi,m = setting$m,
                        maxit = setting$maxit, imputation_method = setting$imputation_method, train_prop = setting$train_prop,
                        val_prop = setting$val_prop,mse_complete = mse_complete,mse_beta_mean = unname(mse_pooled["mean"]),
                        mse_beta_geo  = unname(mse_pooled["geo"]),mse_beta_svd  = unname(mse_pooled["svd"]),
                        R_beta_mean = unname(R["mean"]), R_beta_geo  = unname(R["geo"]), R_beta_svd  = unname(R["svd"]),beta_within_imp_var = beta_within_imp_var,
                        l2_delta_beta_mean = l2_delta_beta_mean,l2_delta_beta_geo  = l2_delta_beta_geo,l2_delta_beta_svd  = l2_delta_beta_svd)
  
  return(list(results = results,complete_obj = complete_obj,split_obj = split_obj,bench_out = bench_out,
              train_data_missing = train_data_missing,imp_out = imp_out,imputed_out = imputed_out,
              pooled_beta = pooled_beta, mse_pooled = mse_pooled,R = R,beta_within_imp_var = beta_within_imp_var,true_beta_full = true_beta_full))
  }