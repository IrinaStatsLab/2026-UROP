
# simulation/R/data_visualization.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

dir_rds <- "output/rds"
dir_fig <- "output/figures"

if (!dir.exists(dir_fig)) dir.create(dir_fig, recursive = TRUE)

# =========================================================
# 0) Helper functions
# =========================================================

read_one <- function(file, beta_type, mechanism, setting, lambda_mode) {
  readRDS(file) %>%
    mutate(
      beta_type   = beta_type,
      mechanism   = mechanism,
      setting     = setting,      # baseline or train_test
      lambda_mode = lambda_mode   # cv / fixed_1 / train_test_cv
    )
}

clean_common <- function(df) {
  df %>%
    mutate(
      rho = factor(rho, levels = sort(unique(rho))),
      pi = factor(pi, levels = sort(unique(pi))),
      beta_type = factor(beta_type, levels = c("sparse", "dense")),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR"))
    )
}

theme_set(
  theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
)

save_plot <- function(p, filename, width = 12, height = 8) {
  ggsave(
    filename = file.path(dir_fig, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

# =========================================================
# 1) Read baseline results (NO train-test)
#    Purpose:
#    - estimate target quantity
#    - study factors affecting estimation
#    - compare pooling strategies
#    - compare CV lambda vs fixed lambda = 1
# =========================================================

baseline_cv_df <- bind_rows(
  read_one(file.path(dir_rds, "replicate_results_sparse_MCAR.rds"),
           "sparse", "MCAR", "baseline", "cv"),
  read_one(file.path(dir_rds, "replicate_results_sparse_MAR.rds"),
           "sparse", "MAR", "baseline", "cv"),
  read_one(file.path(dir_rds, "replicate_results_dense_MCAR.rds"),
           "dense", "MCAR", "baseline", "cv"),
  read_one(file.path(dir_rds, "replicate_results_dense_MAR.rds"),
           "dense", "MAR", "baseline", "cv")
) %>%
  clean_common()

baseline_fixed_df <- bind_rows(
  read_one(file.path(dir_rds, "replicate_results_sparse_MCAR_lambda_1.0.rds"),
           "sparse", "MCAR", "baseline", "fixed_1"),
  read_one(file.path(dir_rds, "replicate_results_sparse_MAR_lambda_1.0.rds"),
           "sparse", "MAR", "baseline", "fixed_1"),
  read_one(file.path(dir_rds, "replicate_results_dense_MCAR_lambda_1.0.rds"),
           "dense", "MCAR", "baseline", "fixed_1"),
  read_one(file.path(dir_rds, "replicate_results_dense_MAR_lambda_1.0.rds"),
           "dense", "MAR", "baseline", "fixed_1")
) %>%
  clean_common()

baseline_all_df <- bind_rows(baseline_cv_df, baseline_fixed_df)

# Long format for pooling strategy comparison
baseline_long <- baseline_all_df %>%
  select(
    n, rho, pi, beta_type, mechanism, lambda_mode,
    ratio_R,
    ratio_R_coef_mean,
    ratio_R_coef_svd
  ) %>%
  pivot_longer(
    cols = c(ratio_R, ratio_R_coef_mean, ratio_R_coef_svd),
    names_to = "method",
    values_to = "R"
  ) %>%
  mutate(
    method = recode(
      method,
      ratio_R = "Pool MSE",
      ratio_R_coef_mean = "Pooled Coefficient (Mean)",
      ratio_R_coef_svd = "Pooled Coefficient (SVD)"
    ),
    lambda_mode = recode(
      lambda_mode,
      cv = "CV-selected lambda",
      fixed_1 = "Fixed lambda = 1"
    )
  )

# =========================================================
# 2) Read train-test results
#    Purpose:
#    - compare out-of-sample predictive performance
# =========================================================

train_test_df <- bind_rows(
  read_one(file.path(dir_rds, "replicate_results_train_test_sparse_MCAR.rds"),
           "sparse", "MCAR", "train_test", "train_test_cv"),
  read_one(file.path(dir_rds, "replicate_results_train_test_sparse_MAR.rds"),
           "sparse", "MAR", "train_test", "train_test_cv"),
  read_one(file.path(dir_rds, "replicate_results_train_test_dense_MCAR.rds"),
           "dense", "MCAR", "train_test", "train_test_cv"),
  read_one(file.path(dir_rds, "replicate_results_train_test_dense_MAR.rds"),
           "dense", "MAR", "train_test", "train_test_cv")
) %>%
  clean_common()

train_test_long <- train_test_df %>%
  select(
    n, rho, pi, beta_type, mechanism,
    ratio_R_test,
    ratio_R_coef_mean_test,
    ratio_R_coef_svd_test
  ) %>%
  pivot_longer(
    cols = c(ratio_R_test, ratio_R_coef_mean_test, ratio_R_coef_svd_test),
    names_to = "method",
    values_to = "R"
  ) %>%
  mutate(
    method = recode(
      method,
      ratio_R_test = "Pool MSE",
      ratio_R_coef_mean_test = "Pooled Coefficient (Mean)",
      ratio_R_coef_svd_test = "Pooled Coefficient (SVD)"
    )
  )

# =========================================================
# Part 1: Baseline estimate-level figures
# =========================================================
# =========================================================
# Figure 0: Histogram of R (baseline, CV lambda)
# Purpose: overall distribution of estimator
# =========================================================

ratios <- baseline_cv_df$ratio_R
ratios <- ratios[is.finite(ratios)]

sd_R <- sd(ratios, na.rm = TRUE)

p0 <- ggplot(data.frame(R = ratios), aes(x = R)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 1, sd = sd_R),
    color = "red",
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Overall Distribution of Relative MSE (R)",
    subtitle = "Baseline setting with CV-selected lambda",
    x = "R",
    y = "Density"
  )

print(p0)
save_plot(p0, "fig0_histogram_R_baseline_cv.png", width = 10, height = 6)

# ---------------------------------------------------------
# Figure 1:
# Main effect of missing proportion (baseline, Pool MSE only)
# This is for: what factors affect estimation?
# ---------------------------------------------------------

fig1_df <- baseline_cv_df %>%
  group_by(n, rho, pi, beta_type, mechanism) %>%
  summarise(
    R_mean = mean(ratio_R, na.rm = TRUE),
    .groups = "drop"
  )

p1 <- ggplot(fig1_df, aes(x = pi, y = R_mean, color = rho, group = rho)) +
  geom_boxplot(
    data = baseline_cv_df,
    aes(x = pi, y = ratio_R, group = interaction(pi, rho), fill = rho),
    inherit.aes = FALSE,
    alpha = 0.20,
    position = position_dodge(width = 0.75),
    outlier.size = 0.6
  ) +
  geom_point(size = 2.2, position = position_dodge(width = 0.75)) +
  geom_line(linewidth = 0.8, aes(group = rho)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + n ~ mechanism) +
  labs(
    title = "Main Effect of Missingness on Relative MSE",
    subtitle = "Baseline setting, Pool MSE, CV-selected lambda",
    x = "Missingness proportion (pi)",
    y = "Mean R",
    color = "Collinearity (rho)",
    fill = "Collinearity (rho)"
  )

print(p1)
save_plot(p1, "fig1_main_effect_missingness_baseline_cv.png", width = 13, height = 9)

# ---------------------------------------------------------
# Figure 2:
# Pooling strategy comparison (baseline, CV-selected lambda)
# Teacher prefers boxplots
# ---------------------------------------------------------

baseline_cv_long <- baseline_long %>%
  filter(lambda_mode == "CV-selected lambda")

p2 <- ggplot(baseline_cv_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Distribution of Relative MSE by Pooling Strategy",
    subtitle = "Baseline setting with CV-selected lambda",
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p2)
save_plot(p2, "fig2_pooling_strategy_baseline_cv_boxplot.png", width = 13, height = 9)

# ---------------------------------------------------------
# Figure 3:
# Same comparison but fixed lambda = 1
# ---------------------------------------------------------

baseline_fixed_long <- baseline_long %>%
  filter(lambda_mode == "Fixed lambda = 1")

p3 <- ggplot(baseline_fixed_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Distribution of Relative MSE by Pooling Strategy",
    subtitle = "Baseline setting with fixed lambda = 1",
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p3)
save_plot(p3, "fig3_pooling_strategy_baseline_fixed1_boxplot.png", width = 13, height = 9)

# ---------------------------------------------------------
# Figure 4:
# Direct lambda comparison
# Compare CV-selected lambda vs fixed lambda = 1
# ---------------------------------------------------------

p4 <- ggplot(baseline_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.6
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism + lambda_mode) +
  labs(
    title = "Sensitivity to Lambda Choice",
    subtitle = "Comparison of CV-selected lambda and fixed lambda = 1",
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p4)
save_plot(p4, "fig4_lambda_sensitivity_boxplot.png", width = 16, height = 9)

# ---------------------------------------------------------
# Figure 5:
# Method gap plot
# Positive = coefficient pooling better than direct Pool MSE
# ---------------------------------------------------------

gap_df <- baseline_cv_df %>%
  mutate(
    gap_mean = ratio_R - ratio_R_coef_mean,
    gap_svd  = ratio_R - ratio_R_coef_svd
  ) %>%
  select(n, rho, pi, beta_type, mechanism, gap_mean, gap_svd) %>%
  pivot_longer(
    cols = c(gap_mean, gap_svd),
    names_to = "comparison",
    values_to = "gap"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      gap_mean = "Pool MSE - Pooled Coef (Mean)",
      gap_svd  = "Pool MSE - Pooled Coef (SVD)"
    )
  ) %>%
  clean_common()

p5 <- ggplot(gap_df, aes(x = pi, y = gap, fill = comparison)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Gap Relative to Direct Pooling of MSE",
    subtitle = "Positive values indicate coefficient pooling gives smaller R",
    x = "Missingness proportion (pi)",
    y = "Gap",
    fill = "Comparison"
  )

print(p5)
save_plot(p5, "fig5_gap_vs_poolmse_boxplot.png", width = 13, height = 9)

# =========================================================
# Part 2: Train-test figures
# =========================================================

# ---------------------------------------------------------
# Figure 6:
# Out-of-sample comparison of pooling strategies
# ---------------------------------------------------------

p6 <- ggplot(train_test_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Held-out Test Relative MSE by Pooling Strategy",
    subtitle = "Train-test split setting",
    x = "Missingness proportion (pi)",
    y = "Test Relative MSE (R)",
    fill = "Method"
  )

print(p6)
save_plot(p6, "fig6_train_test_pooling_strategy_boxplot.png", width = 13, height = 9)

# ---------------------------------------------------------
# Figure 7:
# Out-of-sample method gap
# ---------------------------------------------------------

train_gap_df <- train_test_df %>%
  mutate(
    gap_mean = ratio_R_test - ratio_R_coef_mean_test,
    gap_svd  = ratio_R_test - ratio_R_coef_svd_test
  ) %>%
  select(n, rho, pi, beta_type, mechanism, gap_mean, gap_svd) %>%
  pivot_longer(
    cols = c(gap_mean, gap_svd),
    names_to = "comparison",
    values_to = "gap"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      gap_mean = "Pool MSE - Pooled Coef (Mean)",
      gap_svd  = "Pool MSE - Pooled Coef (SVD)"
    )
  ) %>%
  clean_common()

p7 <- ggplot(train_gap_df, aes(x = pi, y = gap, fill = comparison)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Held-out Test Gap Relative to Direct Pooling of MSE",
    subtitle = "Positive values indicate coefficient pooling gives smaller test R",
    x = "Missingness proportion (pi)",
    y = "Gap on test set",
    fill = "Comparison"
  )

print(p7)
save_plot(p7, "fig7_train_test_gap_boxplot.png", width = 13, height = 9)

# =========================================================
# Optional appendix figure:
# Variance / stability
# =========================================================

var_plot_df <- baseline_cv_df %>%
  select(
    n, rho, pi, beta_type, mechanism,
    ratio_R_var,
    ratio_R_coef_mean_var,
    ratio_R_coef_svd_var
  ) %>%
  pivot_longer(
    cols = c(ratio_R_var, ratio_R_coef_mean_var, ratio_R_coef_svd_var),
    names_to = "method",
    values_to = "R_var"
  ) %>%
  mutate(
    method = recode(
      method,
      ratio_R_var = "Pool MSE",
      ratio_R_coef_mean_var = "Pooled Coefficient (Mean)",
      ratio_R_coef_svd_var = "Pooled Coefficient (SVD)"
    )
  ) %>%
  clean_common()

p8 <- ggplot(var_plot_df, aes(x = pi, y = R_var, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.size = 0.7
  ) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Stability of Relative MSE by Pooling Strategy",
    subtitle = "Baseline setting with CV-selected lambda",
    x = "Missingness proportion (pi)",
    y = "Variance of R",
    fill = "Method"
  )

print(p8)
save_plot(p8, "fig8_stability_variance_boxplot.png", width = 13, height = 9)

# =========================================================
# Console summary
# =========================================================

cat("\nAll figures saved to:", dir_fig, "\n")
cat("Main recommended figures:\n")
cat("  fig0_histogram_R_baseline_cv.png\n")
cat("  fig1_main_effect_missingness_baseline_cv.png\n")
cat("  fig2_pooling_strategy_baseline_cv_boxplot.png\n")
cat("  fig4_lambda_sensitivity_boxplot.png\n")
cat("  fig6_train_test_pooling_strategy_boxplot.png\n\n")