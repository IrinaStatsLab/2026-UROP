
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
      setting     = setting,
      lambda_mode = lambda_mode
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
      strip.text = element_text(face = "bold")
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
# 1) Read baseline results
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
# Part 1: Baseline estimation-level figures
# =========================================================

# Figure 0a: Histogram of R for Pool MSE with CV lambda

ratios_cv <- baseline_cv_df$ratio_R
ratios_cv <- ratios_cv[is.finite(ratios_cv)]
sd_R_cv <- sd(ratios_cv, na.rm = TRUE)

p0a <- ggplot(data.frame(R = ratios_cv), aes(x = R)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 1, sd = sd_R_cv),
    color = "red",
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "R",
    y = "Density"
  )

print(p0a)
save_plot(p0a, "fig0a_histogram_R_poolMSE_cv.png", width = 10, height = 6)

# Figure 0b: Histogram of R for Pool MSE with fixed lambda

ratios_fixed <- baseline_fixed_df$ratio_R
ratios_fixed <- ratios_fixed[is.finite(ratios_fixed)]
sd_R_fixed <- sd(ratios_fixed, na.rm = TRUE)

p0b <- ggplot(data.frame(R = ratios_fixed), aes(x = R)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 1, sd = sd_R_fixed),
    color = "red",
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "R",
    y = "Density"
  )

print(p0b)
save_plot(p0b, "fig0b_histogram_R_poolMSE_fixed1.png", width = 10, height = 6)

# Figure 0c: Histogram of R for coefficient pooling (mean)

ratios_coef_mean <- baseline_cv_df$ratio_R_coef_mean
ratios_coef_mean <- ratios_coef_mean[is.finite(ratios_coef_mean)]
sd_R_coef_mean <- sd(ratios_coef_mean, na.rm = TRUE)

p0c <- ggplot(data.frame(R = ratios_coef_mean), aes(x = R)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 1, sd = sd_R_coef_mean),
    color = "red",
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "R",
    y = "Density"
  )

print(p0c)
save_plot(p0c, "fig0c_histogram_R_coefMean_cv.png", width = 10, height = 6)

# Figure 1: Main effect of missingness on Pool MSE

p1 <- ggplot(
  baseline_cv_df,
  aes(x = pi, y = ratio_R)
) +
  geom_boxplot(
    fill = "steelblue",
    color = "black",
    alpha = 0.8,
    outlier.size = 0.6
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + n ~ mechanism + rho) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)"
  )

print(p1)
save_plot(p1, "fig1_main_effect_missingness_baseline_cv.png", width = 15, height = 9)

# Figure 2: Pooling strategy comparison under CV lambda

baseline_cv_long <- baseline_long %>%
  filter(lambda_mode == "CV-selected lambda")

p2 <- ggplot(baseline_cv_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.9),
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 1.3)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p2)
save_plot(p2, "fig2_pooling_strategy_baseline_cv_boxplot.png", width = 13, height = 9)

# Figure 3: Pooling strategy comparison under fixed lambda

baseline_fixed_long <- baseline_long %>%
  filter(lambda_mode == "Fixed lambda = 1")

p3 <- ggplot(baseline_fixed_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0.6, 1.3)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p3)
save_plot(p3, "fig3_pooling_strategy_baseline_fixed1_boxplot.png", width = 13, height = 9)

# Figure 4: Sensitivity to lambda choice

p4 <- ggplot(baseline_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism + lambda_mode) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Relative MSE (R)",
    fill = "Method"
  )

print(p4)
save_plot(p4, "fig4_lambda_sensitivity_boxplot.png", width = 16, height = 9)

# Figure 5: Gap analysis (baseline)

gap_df <- baseline_cv_df %>%
  mutate(
    gap_poolmse_vs_mean = abs(ratio_R - 1) - abs(ratio_R_coef_mean - 1),
    gap_poolmse_vs_svd  = abs(ratio_R - 1) - abs(ratio_R_coef_svd - 1),
    gap_mean_vs_svd     = abs(ratio_R_coef_mean - 1) - abs(ratio_R_coef_svd - 1)
  ) %>%
  clean_common()

p5a <- ggplot(gap_df, aes(x = pi, y = gap_poolmse_vs_mean)) +
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p5a)
save_plot(p5a, "fig5a_gap_poolmse_vs_mean.png", width = 13, height = 9)

p5b <- ggplot(gap_df, aes(x = pi, y = gap_poolmse_vs_svd)) +
  geom_boxplot(
    fill = "darkorange",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p5b)
save_plot(p5b, "fig5b_gap_poolmse_vs_svd.png", width = 13, height = 9)

p5c <- ggplot(gap_df, aes(x = pi, y = gap_mean_vs_svd)) +
  geom_boxplot(
    fill = "mediumpurple",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-0.03, 0.01)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p5c)
save_plot(p5c, "fig5c_gap_mean_vs_svd.png", width = 13, height = 9)

# =========================================================
# Part 2: Train-test figures
# =========================================================

# Figure 6: Out-of-sample comparison

p6 <- ggplot(train_test_long, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(
    alpha = 0.75,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(1, 6)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Test Relative MSE (R)",
    fill = "Method"
  )

print(p6)
save_plot(p6, "fig6_train_test_pooling_strategy_boxplot.png", width = 13, height = 9)

# Figure 6b: Histogram of selected lambda across imputed datasets

lambda_imp_vals <- train_test_df$lambda_imp_mean
lambda_imp_vals <- lambda_imp_vals[is.finite(lambda_imp_vals)]

p6b <- ggplot(data.frame(lambda = lambda_imp_vals), aes(x = lambda)) +
  geom_histogram(
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  labs(
    x = "Selected lambda",
    y = "Count"
  )

print(p6b)
save_plot(p6b, "fig6b_train_test_lambda_histogram.png", width = 10, height = 6)

# Figure 6c: Overlaid histogram of benchmark vs imputed-data lambda

lambda_compare_hist_df <- train_test_df %>%
  select(lambda_full, lambda_imp_mean) %>%
  pivot_longer(
    cols = c(lambda_full, lambda_imp_mean),
    names_to = "lambda_type",
    values_to = "lambda"
  ) %>%
  mutate(
    lambda_type = recode(
      lambda_type,
      lambda_full = "Complete-data benchmark",
      lambda_imp_mean = "Mean across imputed datasets"
    )
  ) %>%
  filter(is.finite(lambda))

p6c <- ggplot(lambda_compare_hist_df, aes(x = lambda, fill = lambda_type)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 30
  ) +
  labs(
    x = "Selected lambda",
    y = "Count",
    fill = "Lambda source"
  )

print(p6c)
save_plot(p6c, "fig6c_train_test_lambda_compare_histogram.png", width = 10, height = 6)

# Figure 7: Gap analysis (train-test)

train_gap_df <- train_test_df %>%
  mutate(
    gap_poolmse_vs_mean = abs(ratio_R_test - 1) - abs(ratio_R_coef_mean_test - 1),
    gap_poolmse_vs_svd  = abs(ratio_R_test - 1) - abs(ratio_R_coef_svd_test - 1),
    gap_mean_vs_svd     = abs(ratio_R_coef_mean_test - 1) - abs(ratio_R_coef_svd_test - 1)
  ) %>%
  clean_common()

p7a <- ggplot(train_gap_df, aes(x = pi, y = gap_poolmse_vs_mean)) +
  geom_boxplot(
    fill = "steelblue",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 4)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p7a)
save_plot(p7a, "fig7a_train_test_gap_poolmse_vs_mean.png", width = 13, height = 9)

p7b <- ggplot(train_gap_df, aes(x = pi, y = gap_poolmse_vs_svd)) +
  geom_boxplot(
    fill = "darkorange",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 3.5)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p7b)
save_plot(p7b, "fig7b_train_test_gap_poolmse_vs_svd.png", width = 13, height = 9)

p7c <- ggplot(train_gap_df, aes(x = pi, y = gap_mean_vs_svd)) +
  geom_boxplot(
    fill = "mediumpurple",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-1, 0)) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    x = "Missingness proportion (pi)",
    y = "Gap (distance to R = 1)"
  )

print(p7c)
save_plot(p7c, "fig7c_train_test_gap_mean_vs_svd.png", width = 13, height = 9)

# =========================================================
# Console summary
# =========================================================

cat("\nAll figures saved to:", dir_fig, "\n")