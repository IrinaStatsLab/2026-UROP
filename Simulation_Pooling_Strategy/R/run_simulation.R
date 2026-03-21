# Simulation_Pooling_Strategy/R/run_simulation.R
# Run full simulation across scenarios and replicates

## Load packages
library(ggplot2)
library(MASS)
library(mice)
library(glmnet)
library(dplyr)
library(foreach)
library(doParallel)
library(doRNG)

## Source helper files
source("Simulation_Pooling_Strategy/R/config_protocol.R")
source("Simulation_Pooling_Strategy/R/make_setting.R")
source("Simulation_Pooling_Strategy/R/data_generation.R")
source("Simulation_Pooling_Strategy/R/split_data.R")
source("Simulation_Pooling_Strategy/R/impose_missingness.R")
source("Simulation_Pooling_Strategy/R/run_mice_imputation.R")
source("Simulation_Pooling_Strategy/R/fit_complete_dataset.R")
source("Simulation_Pooling_Strategy/R/fit_imputed_datasets.R")
source("Simulation_Pooling_Strategy/R/pooling_coefficients.R")
source("Simulation_Pooling_Strategy/R/evaluate_results.R")
source("Simulation_Pooling_Strategy/R/one_simulation_replicate.R")

## Run simulation
run_simulation <- function(n_grid, rho_grid, pi_grid, beta_grid, mech_grid,
                           R_reps, seed0, ncores,save_plots = TRUE, 
                           plot_dir = "Simulation_Pooling_Strategy/outputs") {
  
  scenarios <- expand.grid(n = n_grid,rho = rho_grid, pi = pi_grid, beta_type = beta_grid, mechanism = mech_grid, stringsAsFactors = FALSE)
  
  ## sanitize ncores
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
  
  ## Store all replicate-level results and scenario summaries
  results_all <- vector("list", nrow(scenarios))
  summary_all <- vector("list", nrow(scenarios))
  
  ## functions / globals workers need
  export_objs <- c("make_setting", "run_one_simulation", "generate_complete_data", "make_sigma_block_ar1",
                   "make_beta", "inv_logit", "split_data", "impose_mcar", "impose_mar", "impose_missingness",
                   "run_mice_imputation", "fit_complete_benchmark", "fit_imputed_datasets", "pool_beta_mean", "pool_beta_geo", "pool_beta_svd", "pool_beta_all",
                   "predict_from_beta", "compute_mse", "evaluate_pooled_betas", "compute_relative_mse")
  
  export_globals <- c("b", "K", "s", "signal", "sigma_eps","train_prop", "val_prop", "M_risk", "a1", "a2", "a3",
                      "m", "maxit", "imputation_method","lambda_grid", "seed0")
  
  for (scen in seq_len(nrow(scenarios))) {
    sc <- scenarios[scen, ]
    
    cat(
      "\nRunning scenario:",
      "n =", sc$n,
      ", rho =", sc$rho,
      ", pi =", sc$pi,
      ", beta_type =", sc$beta_type,
      ", mechanism =", sc$mechanism, "\n"
    )
    
    ## Run Monte Carlo replications for this scenario in parallel
    res <- foreach::foreach(
      r = seq_len(R_reps),
      .combine = dplyr::bind_rows,
      .packages = c("MASS", "mice", "glmnet", "dplyr"),
      .export = c(export_objs, export_globals)
    ) %dopar% {
      setting <- make_setting(
        n = sc$n,
        rho = sc$rho,
        pi = sc$pi,
        beta_type = sc$beta_type,
        mechanism = sc$mechanism,
        rep_id = r
      )
      run_one_simulation(setting)$results
    }
    
    ## Store replicate-level results
    results_all[[scen]] <- res
    
    ## Scenario summary
    summ <- res %>%
      summarise(
        n = first(n),
        rho = first(rho),
        pi = first(pi),
        beta_type = first(beta_type),
        mechanism = first(mechanism),
        mse_complete_mean = mean(mse_complete, na.rm = TRUE),
        mse_beta_mean_mean = mean(mse_beta_mean, na.rm = TRUE),
        mse_beta_geo_mean  = mean(mse_beta_geo, na.rm = TRUE),
        mse_beta_svd_mean  = mean(mse_beta_svd, na.rm = TRUE),
        R_beta_mean_mean = mean(R_beta_mean, na.rm = TRUE),
        R_beta_mean_sd   = sd(R_beta_mean, na.rm = TRUE),
        R_beta_geo_mean  = mean(R_beta_geo, na.rm = TRUE),
        R_beta_geo_sd    = sd(R_beta_geo, na.rm = TRUE),
        R_beta_svd_mean  = mean(R_beta_svd, na.rm = TRUE),
        R_beta_svd_sd    = sd(R_beta_svd, na.rm = TRUE),
        beta_within_imp_var_mean = mean(beta_within_imp_var, na.rm = TRUE),
        l2_delta_beta_mean_mean = mean(l2_delta_beta_mean, na.rm = TRUE),
        l2_delta_beta_geo_mean  = mean(l2_delta_beta_geo, na.rm = TRUE),
        l2_delta_beta_svd_mean  = mean(l2_delta_beta_svd, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(summ)
    summary_all[[scen]] <- summ
  }
  
  ## Combine all outputs
  final_results <- dplyr::bind_rows(results_all)
  summary_results <- dplyr::bind_rows(summary_all)
  
  ## Save outputs
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  write.csv(final_results, file = file.path(plot_dir, "final_results_all_reps.csv"), row.names = FALSE)
  
  write.csv(summary_results,file = file.path(plot_dir, "summary_results.csv"),row.names = FALSE)
  
  return(list(final_results = final_results, summary_results = summary_results))}

## Execute
cat("\n--- Running pooling-strategy simulation ---\n")

sim_out <- run_simulation(n_grid = n_grid,rho_grid = rho_grid,pi_grid = pi_grid, beta_grid = beta_grid,
                          mech_grid = mech_grid,R_reps = R_reps,seed0 = seed0,ncores = ncores)
results <- sim_out$summary_results

## Data visualization
res_plot <- results %>%
  mutate(rho = factor(rho, levels = sort(unique(rho))),pi  = factor(pi, levels = sort(unique(pi))),
         beta_type = factor(beta_type, levels = c("sparse", "dense")),mechanism = factor(mechanism, levels = c("MCAR", "MAR")))

## Heatmap for R_beta_mean
p_heat_mean <- ggplot(res_plot, aes(x = pi, y = rho, fill = R_beta_mean_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", R_beta_mean_mean)), size = 3) +
  facet_grid(beta_type ~ mechanism) +
  labs(
    title = "Mean relative MSE for mean pooling",
    x = "Missingness level (pi)",
    y = "Collinearity (rho)",
    fill = "R mean"
  ) +
  theme_minimal()
print(p_heat_mean)

## Heatmap for R_beta_geo
p_heat_geo <- ggplot(res_plot, aes(x = pi, y = rho, fill = R_beta_geo_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", R_beta_geo_mean)), size = 3) +
  facet_grid(beta_type ~ mechanism) +
  labs(
    title = "Mean relative MSE for geometric pooling",
    x = "Missingness level (pi)",
    y = "Collinearity (rho)",
    fill = "R geo"
  ) +
  theme_minimal()
print(p_heat_geo)

## Heatmap for R_beta_svd
p_heat_svd <- ggplot(res_plot, aes(x = pi, y = rho, fill = R_beta_svd_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", R_beta_svd_mean)), size = 3) +
  facet_grid(beta_type ~ mechanism) +
  labs(
    title = "Mean relative MSE for SVD pooling",
    x = "Missingness level (pi)",
    y = "Collinearity (rho)",
    fill = "R svd"
  ) +
  theme_minimal()
print(p_heat_svd)

cat("\n--- DONE. Scenario summaries: ---\n")
results_table <- results %>%
  arrange(mechanism, beta_type, n, rho, pi)
print(results_table)
