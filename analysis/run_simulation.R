
# simulation/R/run_simulation.R
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
source("simulation/R/config_protocol.R") 
source("simulation/R/data_generation.R")
source("simulation/R/impose_missingness.R")
source("simulation/R/analysis_utility.R")
source("simulation/R/pooling_coeffcients.R")
source("simulation/R/one_simulation_replicate.R")

## Run simulation 
run_simulation <- function(rho_grid, pi_grid, beta_grid, mech_grid,
                           R_reps, seed0, ncores, save_plots = TRUE, plot_dir = "simulation/outputs") {
  
  scenarios <- expand.grid(
    n = n_grid,
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
  results_all <- vector("list", nrow(scenarios))
  
  ## objects/functions that each worker must see
  export_objs <- c(
    "run_one_replicate",
    "generate_complete_data",
    "impose_mcar", "impose_mar",
    "fit_ridge_fixed_lambda", "predict_ridge", "coef_ridge",
    "intercept_ridge", "predict_from_pooled_beta",
    "pool_beta_mean",
    "make_sigma_block_ar1", "make_beta",
    "mse", "ridge_select_perm_threshold",
    "run_mice_pmm",
    "inv_logit"
  )
  
  ## globals from config_protocol.R used inside run_one_replicate call
  export_globals <- c(
    "b", "K", "s", "signal", "sigma_eps",
    "M_risk", "lambda_fixed", "m_imp", "perm_m", "sel_alpha"
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
        n = sc$n, b = b, K = K, s = s, signal = signal, sigma_eps = sigma_eps,
        M_risk = M_risk, lambda_fixed = lambda_fixed, m_imp = m_imp,
        perm_m = perm_m, sel_alpha = sel_alpha,
        seed = seed0 + 100000 * scen + r
      )
    }
    
    # MC simulation
    ratios <- res$ratio_R
    mu_hat <- mean(ratios, na.rm = TRUE)
    sd_hat <- sd(ratios, na.rm = TRUE)
    Z_vals <- (ratios - mu_hat) / sd_hat
    Z_vals <- Z_vals[is.finite(Z_vals)]  # drop NA/Inf if sd_hat==0 etc.
    tag_mc <- paste0("n", sc$n, "_rho", sc$rho, "_pi", sc$pi, "_", sc$beta_type, "_", sc$mechanism)
    tag_mc <- gsub("[^A-Za-z0-9_.-]", "_", tag_mc)
    
    if (save_plots) {
      if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
      
      png(file.path(plot_dir, paste0("MC_R_", tag_mc, ".png")),
          width = 1000, height = 700)
      
      hist(ratios, breaks = 40, freq = FALSE,
           main = paste0("Raw Ratio R | ", tag_mc),
           xlab = "Ratio R")
      
      # optional: overlay a normal curve using sample mean/sd
      curve(dnorm(x, mean(ratios, na.rm = TRUE), 
                  sd(ratios, na.rm = TRUE)),
            add = TRUE, col = "red", lwd = 2)
      
      dev.off()
    }
    
    # MC variance diagnostics
    mse_imp_vals <- res$pooled_mse
    var_imp <- var(mse_imp_vals, na.rm = TRUE)
    delta_var_ratio <- var_imp / (length(na.omit(mse_imp_vals)) * mean(res$mse_full, na.rm=TRUE)^2)
    delta_sd_ratio  <- sqrt(delta_var_ratio)
    
    ## Selection stability plot (Times Selected) per scenario
    sel_cols <- grep("^sel_X\\d+$", names(res), value = TRUE)
    if (save_plots && length(sel_cols) > 0) {
      # Coerce to numeric 0/1 and sum across replicates; NA treated as 0
      sel_counts <- colSums(as.data.frame(lapply(res[sel_cols], function(x) {
        x <- as.numeric(x)
        x[is.na(x)] <- 0
        x
      })))
      
      df_sel <- data.frame(
        Feature = sub("^sel_", "", names(sel_counts)),
        TimesSelected = as.integer(sel_counts),
        stringsAsFactors = FALSE
      )
      df_sel <- df_sel[order(df_sel$TimesSelected, decreasing = TRUE), ]
      
      p_sel <- ggplot(df_sel, aes(x = TimesSelected, y = reorder(Feature, TimesSelected))) +
        geom_col() +
        labs(
          title = paste0(
            "Feature Selection Stability (", R_reps, " runs) | ",
            "rho=", sc$rho, ", pi=", sc$pi, ", ", sc$beta_type, ", ", sc$mechanism
          ),
          x = "Times Selected",
          y = "Feature"
        ) +
        theme_minimal()
      
      tag <- paste0("n", sc$n, "_rho", sc$rho, "_pi", sc$pi, "_", sc$beta_type, "_", sc$mechanism)
      tag <- gsub("[^A-Za-z0-9_.-]", "_", tag)
      ggplot2::ggsave(
        filename = file.path(plot_dir, paste0("stability_", tag, ".png")),
        plot = p_sel, width = 11, height = 7, dpi = 200
      )
      utils::write.csv(
        df_sel,
        file = file.path(plot_dir, paste0("stability_", tag, ".csv")),
        row.names = FALSE
      )
    }
    
    summ <- res %>%
      summarise(
        n = first(n),
        rho = first(rho),
        pi = first(pi),
        beta_type = first(beta_type),
        mechanism = first(mechanism),
        mse_full_mean = mean(mse_full, na.rm = TRUE),
        
        pooled_mse_mean = mean(pooled_mse, na.rm = TRUE),
        pooled_mse_sd   = sd(pooled_mse, na.rm = TRUE),
        
        ratio_R_mean    = mean(ratio_R, na.rm = TRUE),
        ratio_R_sd      = sd(ratio_R, na.rm = TRUE),
        
        pooled_coef_mse_mean = mean(pooled_coef_mse, na.rm = TRUE),
        pooled_coef_mse_sd   = sd(pooled_coef_mse, na.rm = TRUE),
        
        ratio_R_coef_mean = mean(ratio_R_coef, na.rm = TRUE),
        ratio_R_coef_sd   = sd(ratio_R_coef, na.rm = TRUE),
        
        mse_diff_mean = mean(mse_diff, na.rm = TRUE),
        mse_diff_sd   = sd(mse_diff, na.rm = TRUE),
        
        sel_threshold_mean = if ("sel_threshold" %in% names(res)) mean(sel_threshold, na.rm = TRUE) else NA_real_,
        n_selected_mean = if ("n_selected" %in% names(res)) mean(n_selected, na.rm = TRUE) else NA_real_,
        imp_fail_rate = mean(imp_failed),
        .groups = "drop"
      )
    
    print(summ)
    
    out_all[[scen]] <- summ
    results_all[[scen]] <- res
  }
  
  ## return scenario summaries 
  list(
    summary = dplyr::bind_rows(out_all),
    replicate = dplyr::bind_rows(results_all)
  )
}

## Execute
cat("\n--- Running paired simulation ---\n")
sim_out <- run_simulation(
  rho_grid = rho_grid,
  pi_grid = pi_grid,
  beta_grid = beta_grid,
  mech_grid = mech_grid,
  R_reps = R_reps,
  seed0 = seed0,
  ncores = ncores,
  save_plots = FALSE
)
results <- sim_out$summary
replicate_results <- sim_out$replicate

library(tidyr)

# Data Visualization: compare Pool MSE vs Pooled Coefficient across missingness level
res_plot <- results %>%
  mutate(
    rho = factor(rho, levels = sort(unique(rho))),
    pi  = factor(pi, levels = sort(unique(pi))),
    beta_type = factor(beta_type, levels = c("sparse", "dense")),
    mechanism = factor(mechanism, levels = c("MCAR", "MAR"))
  ) %>%
  pivot_longer(
    cols = c(ratio_R_mean, ratio_R_coef_mean),
    names_to = "method",
    values_to = "R_mean"
  ) %>%
  mutate(
    method = recode(
      method,
      ratio_R_mean = "Pool MSE",
      ratio_R_coef_mean = "Pooled Coefficient"
    )
  )

p_compare <- ggplot(
  res_plot,
  aes(x = pi, y = R_mean, group = method, color = method, shape = method)
) +
  geom_point(size = 2.8) +
  geom_line(linewidth = 0.8) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Relative MSE Compared with Complete-Data Benchmark",
    x = "Missingness proportion (pi)",
    y = "R = MSE on imputed data / MSE on complete data",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 14)

print(p_compare)

cat("\n--- DONE. Scenario summaries: ---\n")
results_table <- results %>%
  arrange(mechanism, beta_type, n, rho, pi)
print(results_table)