
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
                              M_risk, lambda_fixed, m_imp, seed) {
   set.seed(seed)
   
   p <- b * K
  
  ## Generate one complete dataset
   dat <- generate_complete_data(
     n = n, b = b, K = K, rho = rho, beta_type = beta_type,
     s = s, signal = signal, sigma_eps = sigma_eps, seed = seed
   )
   X_full <- dat$X
   y <- dat$y
   
  ## Fit ridge model on complete dataset
  fit_full <- fit_ridge_fixed_lambda(X_full, y, lambda_fixed)
  yhat_full <- predict_ridge(fit_full, X_full)
  mse_full <- mse(y, yhat_full)
  # Store oracle coefficients for diagnostic comparison
  beta_full <- as.numeric(coef(fit_full))[-1]
  
  ## Impose missingness on a COPY of X_full 
  X_mis <- if (mechanism == "MCAR") {
    impose_mcar(X_full, pi, M_risk, seed = seed + 10000)
  } else {
    impose_mar(X_full, y, pi, M_risk, seed = seed + 20000)
  }
  
  ## Multiple imputation on the incomplete predictors
  df_mis <- as.data.frame(X_mis)
  df_mis$y <- y
  
  # Use the wrapper so a single failure doesn't crash the simulation
  imp <- run_mice_pmm(df_mis, m = m_imp, seed = seed + 30000)
  
  # If MICE fails, return NA summaries for missing arm 
  if (is.null(imp)) {
    return(data.frame(
      rho=rho, pi=pi, beta_type=beta_type, mechanism=mechanism,
      mse_full=mse_full, pooled_mse=NA_real_, ratio_R=NA_real_,
      mse_imp_sd=NA_real_, within_imp_var=NA_real_,
      sel_stability_jaccard=NA_real_, l2_delta_beta=NA_real_,
      imp_failed=TRUE
    ))
  }
  
  # Containers for per-imputation results
  mse_imp <- rep(NA_real_, m_imp) # MSE per imputation
  sel_list <- vector("list", m_imp) # top-k index set per imputation
  beta_imp_mat <- matrix(NA_real_, nrow = p, ncol = m_imp) # beta_hat per imputation
  
  # For each imputed dataset:
    # 1. fit ridge
    # 2. evaluate on SAME complete X_full 
    # 3. store MSE and coefficient diagnostics
  
  for (j in seq_len(m_imp)) {
    df_j <- tryCatch(complete(imp, action = j), error = function(e) NULL)
    if (is.null(df_j)) next
    
    ## Extract predictors 
    Xj <- as.matrix(df_j[, paste0("X", seq_len(p)), drop = FALSE])
    fit_j <- fit_ridge_fixed_lambda(Xj, y, lambda_fixed)
    
    ## Evaluate on SAME complete data 
    yhat_j <- predict_ridge(fit_j, X_full)
    mse_imp[j] <- mse(y, yhat_j)
    
    beta_hat <- as.numeric(coef(fit_j))[-1]
    beta_imp_mat[, j] <- beta_hat
    # Ridge doesn't "select" variables, but top-k acts as a stability proxy
    sel_list[[j]] <- selected_topk(beta_hat, s_select = s)
  }
  
  # Pool Across Imputations
  pooled_mse <- mean(mse_imp, na.rm = TRUE)
  ratio_R <- pooled_mse / mse_full
  
  ## Variability of MSE across imputations 
  mse_imp_sd <- sd(mse_imp, na.rm = TRUE)
  within_imp_var <- var(mse_imp, na.rm = TRUE)
  
  #  Selection stability diagnostic 
  jac_vals <- c()
  good_idx <- which(!sapply(sel_list, is.null))
  if (length(good_idx) >= 2) {
    for (a in 1:(length(good_idx) - 1)) for (b2 in (a + 1):length(good_idx)) {
      jac_vals <- c(jac_vals, jaccard(sel_list[[good_idx[a]]], sel_list[[good_idx[b2]]]))
    }
  }
  sel_stability <- mean(jac_vals, na.rm = TRUE)
  
  # Coefficient diagnostic
  beta_imp_mean <- rowMeans(beta_imp_mat, na.rm = TRUE)
  delta_beta <- beta_imp_mean - beta_full
  l2_delta_beta <- sqrt(sum(delta_beta^2))
  
  # Return one-row summary for this replicate
  data.frame(
    rho=rho, pi=pi, beta_type=beta_type, mechanism=mechanism,
    mse_full=mse_full, pooled_mse=pooled_mse, ratio_R=ratio_R,
    mse_imp_sd=mse_imp_sd, within_imp_var=within_imp_var,
    sel_stability_jaccard=sel_stability, l2_delta_beta=l2_delta_beta,
    imp_failed=FALSE
  )
}
