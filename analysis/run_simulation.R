
# simulation/R/run_simulation.R
# Run full simulation across scenarios and replicates

library(MASS)
library(mice)
library(glmnet)
library(dplyr)
library(foreach)
library(doParallel)
library(doRNG)

source("R/simulation_setting.R")
source("R/data_generation.R")
source("R/impose_missingness.R")
source("R/analysis_utility.R")
source("R/pooling_coeffcients.R")
source("R/one_simulation_replicate.R")

run_simulation <- function(rho_grid, pi_grid, beta_grid, mech_grid,
                           R_reps, seed0, ncores) {
  
  scenarios <- expand.grid(
    n = n_grid,
    rho = rho_grid,
    pi = pi_grid,
    beta_type = beta_grid,
    mechanism = mech_grid,
    stringsAsFactors = FALSE
  )
  
  ncores <- as.integer(ncores)
  if (is.na(ncores) || ncores < 1) ncores <- 1
  max_cores <- parallel::detectCores()
  ncores <- min(ncores, max_cores)
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed0)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  results_all <- vector("list", nrow(scenarios))
  
  export_objs <- c(
    "run_one_replicate",
    "generate_complete_data",
    "impose_mcar", "impose_mar",
    "fit_ridge_fixed_lambda", "predict_ridge", "coef_ridge",
    "intercept_ridge", "predict_from_pooled_beta",
    "pool_beta_mean", "beta_pool_svd",
    "make_sigma_block_ar1", "make_beta",
    "mse", "ridge_select_perm_threshold",
    "run_mice_pmm", "select_ridge_lambda_cv",
    "inv_logit"
  )
  
  export_globals <- c(
    "b", "K", "s", "signal", "sigma_eps",
    "M_risk", "lambda_mode", "lambda_fixed_value",
    "m_imp", "perm_m", "sel_alpha"
  )
  
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
    
    results_all[[scen]] <- res
  }
  
  replicate_results <- dplyr::bind_rows(results_all)
  
  results_table <- replicate_results %>%
    group_by(n, rho, pi, beta_type, mechanism) %>%
    summarise(
      lambda_mean = mean(lambda, na.rm = TRUE),
      lambda_sd   = sd(lambda, na.rm = TRUE),
      
      mse_full_mean = mean(mse_full, na.rm = TRUE),
      mse_full_sd   = sd(mse_full, na.rm = TRUE),
      
      pooled_mse_mean = mean(pooled_mse, na.rm = TRUE),
      pooled_mse_sd   = sd(pooled_mse, na.rm = TRUE),
      
      ratio_R_mean = mean(ratio_R, na.rm = TRUE),
      ratio_R_sd   = sd(ratio_R, na.rm = TRUE),
      ratio_R_var_mean = mean(ratio_R_var, na.rm = TRUE),
      
      pooled_coef_mse_mean_mean = mean(pooled_coef_mse_mean, na.rm = TRUE),
      pooled_coef_mse_mean_sd   = sd(pooled_coef_mse_mean, na.rm = TRUE),
      
      ratio_R_coef_mean_mean = mean(ratio_R_coef_mean, na.rm = TRUE),
      ratio_R_coef_mean_sd   = sd(ratio_R_coef_mean, na.rm = TRUE),
      ratio_R_coef_mean_var_mean = mean(ratio_R_coef_mean_var, na.rm = TRUE),
      
      pooled_coef_mse_svd_mean = mean(pooled_coef_mse_svd, na.rm = TRUE),
      pooled_coef_mse_svd_sd   = sd(pooled_coef_mse_svd, na.rm = TRUE),
      
      ratio_R_coef_svd_mean = mean(ratio_R_coef_svd, na.rm = TRUE),
      ratio_R_coef_svd_sd   = sd(ratio_R_coef_svd, na.rm = TRUE),
      ratio_R_coef_svd_var_mean = mean(ratio_R_coef_svd_var, na.rm = TRUE),
      
      sel_threshold_mean = mean(sel_threshold, na.rm = TRUE),
      n_selected_mean    = mean(n_selected, na.rm = TRUE),
      
      l2_delta_beta_mean_mean = mean(l2_delta_beta_mean, na.rm = TRUE),
      l2_delta_beta_svd_mean  = mean(l2_delta_beta_svd, na.rm = TRUE),
      
      imp_fail_rate = mean(imp_failed, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    summary = results_table,
    replicate = replicate_results
  )
}

cat("\n--- Running paired simulation ---\n")

sim_out <- run_simulation(
  rho_grid = rho_grid,
  pi_grid = pi_grid,
  beta_grid = beta_grid,
  mech_grid = mech_grid,
  R_reps = R_reps,
  seed0 = seed0,
  ncores = ncores
)

results_table <- sim_out$summary
replicate_results <- sim_out$replicate

if (!dir.exists("output")) dir.create("output", recursive = TRUE)
if (!dir.exists("output/rds")) dir.create("output/rds", recursive = TRUE)

write.csv(results_table, "output/tables/results_table.csv", row.names = FALSE)
saveRDS(replicate_results, "output/rds/replicate_results.rds")

cat("\n--- DONE. Files saved: ---\n")
cat("Summary table: output/tables/results_table.csv\n")
cat("Replicate-level data: output/rds/replicate_results.rds\n")