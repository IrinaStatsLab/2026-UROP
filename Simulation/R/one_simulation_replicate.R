
# simulation/R/one_simulation_replicate.R
# Purpose: Compare a ridge fit under complete dataset and under the same dataset but with missingness
# Design:
# 1) Generate a single complete dataset (X_full, y).
# 2) Complete Arm: fit ridge on X_full and evaluate on X_full.
# 3) Missing Arm: impose missingness on a COPY of X_full (X_mis), perform multiple imputation via MICE (PMM),
# fit ridge on each imputed dataset, evaluate each fitted model on the SAME X_full, and pool performance by averaging MSE across imputations.

# Outputs: Returns a one-row summary including oracle MSE, pooled MI MSE, their ratio R, and optional diagnostics (MSE variability, stability proxy, etc.).

run_one_replicate <- function(rho, pi, beta_type, mechanism,
                              n, b, K, s, signal, sigma_eps,
                              M_risk, lambda_fixed, m_imp, perm_m, sel_alpha, seed) {
  set.seed(seed)
  p <- b * K
  
  ## Generate one complete dataset
  dat <- generate_complete_data(
    n = n, b = b, K = K, rho = rho, beta_type = beta_type,
    s = s, signal = signal, sigma_eps = sigma_eps, seed = seed
  )
  X_full <- as.matrix(dat$X)
  y <- as.numeric(dat$y)
  colnames(X_full) <- paste0("X", seq_len(p))
  
  ## Complete-data benchmark
  fit_full <- fit_ridge_fixed_lambda(X_full, y, lambda_fixed)
  yhat_full <- predict_ridge(fit_full, X_full)
  mse_full <- mse(y, yhat_full)
  beta_full <- coef_ridge(fit_full)
  
  ## Impose missingness
  X_mis <- if (mechanism == "MCAR") {
    impose_mcar(X_full, pi, M_risk, seed = seed + 10000)
  } else {
    impose_mar(X_full, y, pi, M_risk, seed = seed + 20000)
  }
  X_mis <- as.matrix(X_mis)
  colnames(X_mis) <- colnames(X_full)
  
  ## Multiple imputation
  df_mis <- as.data.frame(X_mis)
  df_mis$y <- y
  imp <- run_mice_pmm(df_mis, m = m_imp, seed = seed + 30000)
  
  empty_sel <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = p))
  colnames(empty_sel) <- paste0("sel_X", seq_len(p))
  
  if (is.null(imp)) {
    out0 <- data.frame(
      n = n, rho = rho, pi = pi, beta_type = beta_type, mechanism = mechanism,
      mse_full = mse_full,
      pooled_mse = NA_real_,              # path 1: pool MSE
      ratio_R = NA_real_,
      pooled_coef_mse = NA_real_,         # path 2: pool coefficient
      ratio_R_coef = NA_real_,
      mse_imp_sd = NA_real_,
      within_imp_var = NA_real_,
      mse_coef_sd = NA_real_,
      within_coef_var = NA_real_,
      sel_threshold = NA_real_,
      n_selected = NA_real_,
      l2_delta_beta = NA_real_,
      imp_failed = TRUE
    )
    return(cbind(out0, empty_sel))
  }
  
  ## Containers
  mse_imp <- rep(NA_real_, m_imp)              # path 1 MSEs
  mse_beta_pool <- rep(NA_real_, m_imp)        # path 2 MSEs
  beta_imp_mat <- matrix(NA_real_, nrow = p, ncol = m_imp)
  intercept_imp <- rep(NA_real_, m_imp)
  sel_mat <- matrix(NA_real_, nrow = p, ncol = m_imp)
  thr_vec <- rep(NA_real_, m_imp)
  X_imp_list <- vector("list", m_imp)
  
  ## Per-imputation fit
  for (j in seq_len(m_imp)) {
    df_j <- tryCatch(mice::complete(imp, action = j), error = function(e) NULL)
    if (is.null(df_j)) next
    
    Xj <- as.matrix(df_j[, paste0("X", seq_len(p)), drop = FALSE])
    colnames(Xj) <- paste0("X", seq_len(p))
    X_imp_list[[j]] <- Xj
    
    fit_j <- fit_ridge_fixed_lambda(Xj, y, lambda_fixed)
    
    ## Path 1: each imputation uses its own coefficients
    yhat_j <- predict_ridge(fit_j, Xj)
    mse_imp[j] <- mse(y, yhat_j)
    
    beta_hat <- coef_ridge(fit_j)
    beta_imp_mat[, j] <- as.numeric(beta_hat[paste0("X", seq_len(p))])
    intercept_imp[j] <- intercept_ridge(fit_j)
    
    sel_res <- ridge_select_perm_threshold(
      X = Xj, y = y, lambda_fixed = lambda_fixed,
      perm_m = perm_m, sel_alpha = sel_alpha,
      seed = seed + 40000 + j
    )
    thr_vec[j] <- sel_res$threshold
    
    sel_vec <- as.integer(sel_res$selected)
    names(sel_vec) <- names(sel_res$selected)
    sel_mat[, j] <- sel_vec[paste0("X", seq_len(p))]
  }
  
  ## Path 2: simple mean pooling for coefficients
  good_beta_cols <- which(colSums(!is.na(beta_imp_mat)) == p)
  if (length(good_beta_cols) > 0) {
    beta_pool_mean <- pool_beta_mean(beta_imp_mat[, good_beta_cols, drop = FALSE])
    intercept_pool_mean <- mean(intercept_imp[good_beta_cols], na.rm = TRUE)
    
    for (j in good_beta_cols) {
      Xj <- X_imp_list[[j]]
      if (is.null(Xj)) next
      yhat_pool_j <- predict_from_pooled_beta(beta_pool_mean, intercept_pool_mean, Xj)
      mse_beta_pool[j] <- mse(y, yhat_pool_j)
    }
  } else {
    beta_pool_mean <- rep(NA_real_, p)
  }
  
  ## Summaries
  pooled_mse <- mean(mse_imp, na.rm = TRUE)               # pool MSE path
  ratio_R <- pooled_mse / mse_full
  
  pooled_coef_mse <- mean(mse_beta_pool, na.rm = TRUE)    # pooled coefficient path
  ratio_R_coef <- pooled_coef_mse / mse_full
  mse_diff <- pooled_coef_mse - pooled_mse
  
  mse_imp_sd <- sd(mse_imp, na.rm = TRUE)
  within_imp_var <- var(mse_imp, na.rm = TRUE)
  
  mse_coef_sd <- sd(mse_beta_pool, na.rm = TRUE)
  within_coef_var <- var(mse_beta_pool, na.rm = TRUE)
  
  delta_beta <- beta_pool_mean - as.numeric(beta_full[paste0("X", seq_len(p))])
  l2_delta_beta <- sqrt(sum(delta_beta^2, na.rm = TRUE))
  
  good_cols <- which(colSums(!is.na(sel_mat)) > 0)
  sel_majority <- rep(NA_real_, p)
  if (length(good_cols) > 0) {
    sel_majority <- as.integer(rowMeans(sel_mat[, good_cols, drop = FALSE], na.rm = TRUE) > 0.5)
  }
  sel_df <- as.data.frame(t(sel_majority))
  colnames(sel_df) <- paste0("sel_X", seq_len(p))
  
  sel_threshold <- mean(thr_vec, na.rm = TRUE)
  n_selected <- sum(sel_majority, na.rm = TRUE)
  
  out <- data.frame(
    n = n, rho = rho, pi = pi, beta_type = beta_type, mechanism = mechanism,
    mse_full = mse_full,
    pooled_mse = pooled_mse,              # path 1
    ratio_R = ratio_R,
    pooled_coef_mse = pooled_coef_mse,    # path 2
    ratio_R_coef = ratio_R_coef,
    mse_diff = mse_diff,
    mse_imp_sd = mse_imp_sd,
    within_imp_var = within_imp_var,
    mse_coef_sd = mse_coef_sd,
    within_coef_var = within_coef_var,
    sel_threshold = sel_threshold,
    n_selected = n_selected,
    l2_delta_beta = l2_delta_beta,
    imp_failed = FALSE
  )
  
  cbind(out, sel_df)
}