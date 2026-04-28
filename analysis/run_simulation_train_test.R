
# R/run_simulation_train_test.R
# Run train-test simulation across scenarios and aggregate results

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
source("R/split_data.R")
source("R/train_test_replicate.R")

run_simulation_train_test <- function(rho_grid, pi_grid, beta_grid, mech_grid,
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
    "run_one_replicate_train_test",
    "split_data",
    "generate_complete_data",
    "impose_mcar", "impose_mar",
    "fit_ridge_fixed_lambda", "predict_ridge", "coef_ridge",
    "intercept_ridge", "predict_from_pooled_beta",
    "pool_beta_mean", "beta_pool_svd",
    "make_sigma_block_ar1", "make_beta",
    "mse", "run_mice_pmm", "select_ridge_lambda_cv",
    "inv_logit"
  )
  
  export_globals <- c(
    "b", "K", "s", "signal", "sigma_eps",
    "M_risk", "m_imp"
  )
  
  # Loop over all simulation scenarios
  for (scen in seq_len(nrow(scenarios))) {
    sc <- scenarios[scen, ]
    
    # Run Monte Carlo replicates in parallel
    res <- foreach::foreach(
      r = seq_len(R_reps),
      .combine = dplyr::bind_rows,
      .packages = c("MASS", "mice", "glmnet", "dplyr", "rsample", "doRNG"),
      .export = c(export_objs, export_globals)
    ) %dopar% {
      run_one_replicate_train_test(
        rho = sc$rho,
        pi = sc$pi,
        beta_type = sc$beta_type,
        mechanism = sc$mechanism,
        n = sc$n,
        b = b, K = K, s = s, signal = signal, sigma_eps = sigma_eps,
        M_risk = M_risk,
        m_imp = m_imp,
        seed = seed0 + 100000 * scen + r
      )
    }
    
    results_all[[scen]] <- res
  }
  
  replicate_results <- dplyr::bind_rows(results_all)
  
  # Aggregate replicate-level results by scenario
  results_table <- replicate_results %>%
    group_by(n, rho, pi, beta_type, mechanism) %>%
    summarise(
      lambda_full_mean = mean(lambda_full, na.rm = TRUE),
      lambda_full_sd   = sd(lambda_full, na.rm = TRUE),
      
      lambda_imp_mean_mean = mean(lambda_imp_mean, na.rm = TRUE),
      lambda_imp_mean_sd   = sd(lambda_imp_mean, na.rm = TRUE),
      lambda_imp_sd_mean   = mean(lambda_imp_sd, na.rm = TRUE),
      
      mse_full_test_mean = mean(mse_full_test, na.rm = TRUE),
      mse_full_test_sd   = sd(mse_full_test, na.rm = TRUE),
      
      pooled_mse_test_mean = mean(pooled_mse_test, na.rm = TRUE),
      pooled_mse_test_sd   = sd(pooled_mse_test, na.rm = TRUE),
      
      ratio_R_test_mean = mean(ratio_R_test, na.rm = TRUE),
      ratio_R_test_sd   = sd(ratio_R_test, na.rm = TRUE),
      
      pooled_coef_mse_mean_test_mean = mean(pooled_coef_mse_mean_test, na.rm = TRUE),
      pooled_coef_mse_mean_test_sd   = sd(pooled_coef_mse_mean_test, na.rm = TRUE),
      
      ratio_R_coef_mean_test_mean = mean(ratio_R_coef_mean_test, na.rm = TRUE),
      ratio_R_coef_mean_test_sd   = sd(ratio_R_coef_mean_test, na.rm = TRUE),
      
      pooled_coef_mse_svd_test_mean = mean(pooled_coef_mse_svd_test, na.rm = TRUE),
      pooled_coef_mse_svd_test_sd   = sd(pooled_coef_mse_svd_test, na.rm = TRUE),
      
      ratio_R_coef_svd_test_mean = mean(ratio_R_coef_svd_test, na.rm = TRUE),
      ratio_R_coef_svd_test_sd   = sd(ratio_R_coef_svd_test, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    summary = results_table,
    replicate = replicate_results
  )
}

cat("\n--- Running train/test simulation ---\n")

sim_out_train_test <- run_simulation_train_test(
  rho_grid = rho_grid,
  pi_grid = pi_grid,
  beta_grid = beta_grid,
  mech_grid = mech_grid,
  R_reps = R_reps,
  seed0 = seed0,
  ncores = ncores
)

# Save results
results_table_train_test <- sim_out_train_test$summary
replicate_results_train_test <- sim_out_train_test$replicate

if (!dir.exists("output")) dir.create("output", recursive = TRUE)
if (!dir.exists("output/rds")) dir.create("output/rds", recursive = TRUE)
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

write.csv(results_table_train_test,
          "output/tables/results_table_train_test.csv",
          row.names = FALSE)

saveRDS(replicate_results_train_test,
        "output/rds/replicate_results_train_test.rds")

cat("\n--- DONE. Files saved: ---\n")
cat("Summary table: output/tables/results_table_train_test.csv\n")
cat("Replicate-level data: output/rds/replicate_results_train_test.rds\n")