
# simulation/R/run_simulation.R
# Run full simulation across scenarios and replicates

## Load packages
library(MASS)
library(mice)
library(glmnet)
library(dplyr)
library(foreach)
library(doParallel)
library(doRNG)

## Source helper files
source("simulation/R/data_generation.R")
source("simulation/R/impose_missingness.R")
source("simulation/R/analysis_utility.R")
source("simulation/R/one_simulation_replicate.R")
source("simulation/R/config_protocol.R") 

## Run simulation 
run_simulation <- function(rho_grid, pi_grid, beta_grid, mech_grid,
                           R_reps, seed0, ncores) {
  
  scenarios <- expand.grid(
    rho = rho_grid,
    pi = pi_grid,
    beta_type = beta_grid,
    mechanism = mech_grid,
    stringsAsFactors = FALSE
  )
  
  ## sanitize ncores
  ncores <- as.integer(ncores)
  if (is.na(ncores) || ncores < 1) ncores <- 1
  max_cores <- parallel::detectCores()
  ncores <- min(ncores, max_cores)
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed0)
  
  ## Ensure cluster is stopped even if error happens
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  ## Store scenario-level summaries
  out_all <- vector("list", nrow(scenarios))
  
  ## objects/functions that each worker must see
  export_objs <- c(
    "run_one_replicate",
    "generate_complete_data",
    "impose_mcar", "impose_mar",
    "fit_ridge_fixed_lambda", "predict_ridge",
    "mse", "selected_topk", "jaccard",
    "run_mice_pmm_safe",
    "inv_logit"
  )
  
  ## globals from config_protocol.R used inside run_one_replicate call
  export_globals <- c(
    "n", "b", "K", "s", "signal", "sigma_eps",
    "M_risk", "lambda_fixed", "m_imp"
  )
  
  for (scen in seq_len(nrow(scenarios))) {
    sc <- scenarios[scen, ]
    
    res <- foreach::foreach(
      r = seq_len(R_reps),
      .combine = dplyr::bind_rows,
      .packages = c("MASS", "mice", "glmnet", "dplyr", "doRNG"),
      .export = c(export_objs, export_globals)
    ) %dopar% {
      run_one_replicate(
        rho = sc$rho,
        pi = sc$pi,
        beta_type = sc$beta_type,
        mechanism = sc$mechanism,
        n = n, b = b, K = K, s = s, signal = signal, sigma_eps = sigma_eps,
        M_risk = M_risk, lambda_fixed = lambda_fixed, m_imp = m_imp,
        seed = seed0 + 100000 * scen + r
      )
    }
    
    summ <- res %>%
      summarise(
        rho = first(rho),
        pi = first(pi),
        beta_type = first(beta_type),
        mechanism = first(mechanism),
        mse_full_mean = mean(mse_full, na.rm = TRUE),
        pooled_mse_mean = mean(pooled_mse, na.rm = TRUE),
        ratio_R_mean = mean(ratio_R, na.rm = TRUE),
        ratio_R_sd = sd(ratio_R, na.rm = TRUE),
        imp_fail_rate = mean(imp_failed),
        .groups = "drop"
      )
    
    print(summ)
    
    out_all[[scen]] <- summ
    # results_all[[scen]] <- res
  }
  
  ## return scenario summaries 
  dplyr::bind_rows(out_all)
}

## Execute
cat("\n--- Running paired simulation ---\n")

results <- run_simulation(
  rho_grid = rho_grid,
  pi_grid = pi_grid,
  beta_grid = beta_grid,
  mech_grid = mech_grid,
  R_reps = R_reps,
  seed0 = seed0,
  ncores = ncores
)

cat("\n--- DONE. Scenario summaries: ---\n")
print(results)

# write.csv(results, "ridge_mi_paired_in_sample_results.csv", row.names = FALSE)