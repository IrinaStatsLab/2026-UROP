# R/train_test_replicate.R

run_one_replicate_train_test <- function(rho, pi, beta_type, mechanism,
                                         n, b, K, s, signal, sigma_eps,
                                         M_risk, m_imp, seed,
                                         train_prop = 0.8) {
  set.seed(seed)
  
  p <- b * K
  
  ## 1) Generate complete data
  dat <- generate_complete_data(
    n = n, b = b, K = K, rho = rho, beta_type = beta_type,
    s = s, signal = signal, sigma_eps = sigma_eps, seed = seed
  )
  
  X_full <- as.matrix(dat$X)
  y_full <- as.numeric(dat$y)
  colnames(X_full) <- paste0("X", seq_len(p))
  
  ## 2) Train/test split
  split_obj <- split_data(dat, train_prop = train_prop, seed = seed + 1000)
  
  train_data <- split_obj$train
  test_data  <- split_obj$test
  
  X_train_full <- as.matrix(train_data[, paste0("X", seq_len(p)), drop = FALSE])
  y_train <- as.numeric(train_data$y)
  
  X_test <- as.matrix(test_data[, paste0("X", seq_len(p)), drop = FALSE])
  y_test <- as.numeric(test_data$y)
  
  ## 3) Complete-data benchmark:
  ##    CV on complete training set, evaluate on complete test set
  lambda_full <- select_ridge_lambda_cv(X_train_full, y_train, seed = seed + 5000)
  fit_full <- fit_ridge_fixed_lambda(X_train_full, y_train, lambda_full)
  
  yhat_full_test <- predict_ridge(fit_full, X_test)
  mse_full_test <- mse(y_test, yhat_full_test)
  
  beta_full <- coef_ridge(fit_full)
  
  ## 4) Impose missingness on TRAIN ONLY
  X_train_mis <- if (mechanism == "MCAR") {
    impose_mcar(X_train_full, pi, M_risk, seed = seed + 10000)
  } else if (mechanism == "MAR") {
    impose_mar(X_train_full, y_train, pi, M_risk, seed = seed + 20000)
  } else {
    stop("Unsupported missingness mechanism.")
  }
  
  X_train_mis <- as.matrix(X_train_mis)
  colnames(X_train_mis) <- paste0("X", seq_len(p))
  
  ## 5) MICE on incomplete training set
  df_train_mis <- as.data.frame(X_train_mis)
  df_train_mis$y <- y_train
  
  imp <- run_mice_pmm(df_train_mis, m = m_imp, seed = seed + 30000)
  
  if (is.null(imp)) {
    return(data.frame(
      n = n, rho = rho, pi = pi,
      beta_type = beta_type, mechanism = mechanism,
      mse_full_test = mse_full_test,
      lambda_full = lambda_full,
      pooled_mse_test = NA_real_,
      ratio_R_test = NA_real_,
      pooled_coef_mse_mean_test = NA_real_,
      ratio_R_coef_mean_test = NA_real_,
      pooled_coef_mse_svd_test = NA_real_,
      ratio_R_coef_svd_test = NA_real_,
      lambda_imp_mean = NA_real_,
      lambda_imp_sd = NA_real_,
      imp_failed = TRUE
    ))
  }
  
  ## 6) Per-imputation fit on imputed training sets
  mse_imp_test <- rep(NA_real_, m_imp)
  lambda_imp <- rep(NA_real_, m_imp)
  intercept_imp <- rep(NA_real_, m_imp)
  beta_imp_mat <- matrix(NA_real_, nrow = p, ncol = m_imp)
  rownames(beta_imp_mat) <- paste0("X", seq_len(p))
  
  for (j in seq_len(m_imp)) {
    df_j <- tryCatch(mice::complete(imp, action = j), error = function(e) NULL)
    if (is.null(df_j)) next
    
    Xj <- as.matrix(df_j[, paste0("X", seq_len(p)), drop = FALSE])
    yj <- as.numeric(df_j$y)
    colnames(Xj) <- paste0("X", seq_len(p))
    if (anyNA(Xj) || anyNA(yj)) next
    
    ## each imputed training set gets its own CV lambda
    lambda_j <- select_ridge_lambda_cv(Xj, yj, seed = seed + 40000 + j)
    lambda_imp[j] <- lambda_j
    
    fit_j <- fit_ridge_fixed_lambda(Xj, yj, lambda_j)
    
    ## test-set MSE for pool-MSE strategy
    yhat_j_test <- predict_ridge(fit_j, X_test)
    mse_imp_test[j] <- mse(y_test, yhat_j_test)
    
    ## store coefficients for coefficient-pooling strategies
    beta_j <- coef_ridge(fit_j)
    beta_imp_mat[, j] <- beta_j[paste0("X", seq_len(p))]
    intercept_imp[j] <- intercept_ridge(fit_j)
  }
  
  ## keep only successful imputations
  good_cols <- which(colSums(!is.na(beta_imp_mat)) == p & is.finite(lambda_imp))
  if (length(good_cols) == 0) {
    return(data.frame(
      n = n, rho = rho, pi = pi,
      beta_type = beta_type, mechanism = mechanism,
      mse_full_test = mse_full_test,
      lambda_full = lambda_full,
      pooled_mse_test = NA_real_,
      ratio_R_test = NA_real_,
      pooled_coef_mse_mean_test = NA_real_,
      ratio_R_coef_mean_test = NA_real_,
      pooled_coef_mse_svd_test = NA_real_,
      ratio_R_coef_svd_test = NA_real_,
      lambda_imp_mean = NA_real_,
      lambda_imp_sd = NA_real_,
      imp_failed = TRUE
    ))
  }
  
  ## 7) Pool MSE across imputations
  pooled_mse_test <- mean(mse_imp_test[good_cols], na.rm = TRUE)
  ratio_R_test <- pooled_mse_test / mse_full_test
  
  ## 8) Pool coefficients
  beta_mat_good <- beta_imp_mat[, good_cols, drop = FALSE]
  
  beta_pool_mean <- pool_beta_mean(beta_mat_good)
  beta_pool_svd_vec <- beta_pool_svd(beta_mat_good)   
  
  intercept_pool <- mean(intercept_imp[good_cols], na.rm = TRUE)
  
  ## evaluate pooled coefficients on SAME complete test set
  yhat_mean_test <- predict_from_pooled_beta(beta_pool_mean, intercept_pool, X_test)
  pooled_coef_mse_mean_test <- mse(y_test, yhat_mean_test)
  ratio_R_coef_mean_test <- pooled_coef_mse_mean_test / mse_full_test
  
  yhat_svd_test <- predict_from_pooled_beta(beta_pool_svd_vec, intercept_pool, X_test)
  pooled_coef_mse_svd_test <- mse(y_test, yhat_svd_test)
  ratio_R_coef_svd_test <- pooled_coef_mse_svd_test / mse_full_test
  
  ## optional coefficient distance from complete-data benchmark
  l2_delta_beta_mean <- sqrt(sum((beta_pool_mean - beta_full[paste0("X", seq_len(p))])^2, na.rm = TRUE))
  l2_delta_beta_svd  <- sqrt(sum((beta_pool_svd_vec - beta_full[paste0("X", seq_len(p))])^2, na.rm = TRUE))
  
  ## 9) Output one-row result
  data.frame(
    n = n, rho = rho, pi = pi,
    beta_type = beta_type, mechanism = mechanism,
    
    mse_full_test = mse_full_test,
    lambda_full = lambda_full,
    
    pooled_mse_test = pooled_mse_test,
    ratio_R_test = ratio_R_test,
    
    pooled_coef_mse_mean_test = pooled_coef_mse_mean_test,
    ratio_R_coef_mean_test = ratio_R_coef_mean_test,
    
    pooled_coef_mse_svd_test = pooled_coef_mse_svd_test,
    ratio_R_coef_svd_test = ratio_R_coef_svd_test,
    
    lambda_imp_mean = mean(lambda_imp[good_cols], na.rm = TRUE),
    lambda_imp_sd = sd(lambda_imp[good_cols], na.rm = TRUE),
    
    l2_delta_beta_mean = l2_delta_beta_mean,
    l2_delta_beta_svd = l2_delta_beta_svd,
    
    imp_failed = FALSE
  )
}