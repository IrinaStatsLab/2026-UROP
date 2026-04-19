
# simulation/R/data_visualization.R

library(dplyr)
library(ggplot2)
library(tidyr)

replicate_results <- readRDS("output/rds/replicate_results.rds")

# Keep valid ratio_R values
ratios <- replicate_results$ratio_R
ratios <- ratios[is.finite(ratios)]

# Empirical variance of ratio_R
var_R <- var(ratios, na.rm = TRUE)
sd_R  <- sqrt(var_R)

# Wider x-range around 1
x_min <- min(1 - 8 * sd_R, min(ratios, na.rm = TRUE))
x_max <- max(1 + 8 * sd_R, max(ratios, na.rm = TRUE))

# Histogram + reference N(1, empirical variance)
p_hist_R <- ggplot(data.frame(ratio_R = ratios), aes(x = ratio_R)) +
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
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
  coord_cartesian(xlim = c(x_min, x_max)) +
  labs(
    title = "Distribution of Ratio R",
    subtitle = "Reference curve: N(1, empirical variance of ratio_R)",
    x = "Ratio R",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

print(p_hist_R)

## 2) Compare three methods using replicate-level means already stored
plot_df <- replicate_results %>%
  mutate(
    rho = factor(rho, levels = sort(unique(rho))),
    pi = factor(pi, levels = sort(unique(pi))),
    beta_type = factor(beta_type, levels = c("sparse", "dense")),
    mechanism = factor(mechanism, levels = c("MCAR", "MAR"))
  ) %>%
  select(
    n, rho, pi, beta_type, mechanism,
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
    )
  )

p_box <- ggplot(plot_df, aes(x = pi, y = R, fill = method)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Distribution of Relative MSE by Method",
    x = "Missingness proportion (pi)",
    y = "R",
    fill = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_box)

## 3) Scenario-level line plot from replicate-level data
summary_plot_df <- plot_df %>%
  group_by(n, rho, pi, beta_type, mechanism, method) %>%
  summarise(
    R_mean = mean(R, na.rm = TRUE),
    .groups = "drop"
  )

p_line <- ggplot(
  summary_plot_df,
  aes(x = pi, y = R_mean, color = method, group = method, shape = method)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Mean Relative MSE by Method",
    x = "Missingness proportion (pi)",
    y = "Mean R",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_line)

## 4) Main effects of factors on Pool MSE (baseline)
factor_df <- replicate_results %>%
  group_by(n, rho, pi, beta_type, mechanism) %>%
  summarise(
    R_mean = mean(ratio_R, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rho = factor(rho, levels = sort(unique(rho))),
    pi = factor(pi, levels = sort(unique(pi))),
    beta_type = factor(beta_type, levels = c("sparse", "dense")),
    mechanism = factor(mechanism, levels = c("MCAR", "MAR"))
  )

p_factor <- ggplot(
  factor_df,
  aes(x = pi, y = R_mean, color = rho, group = rho)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type + n ~ mechanism) +
  labs(
    title = "Main Effects on Relative MSE (Pool MSE Baseline)",
    x = "Missingness proportion (pi)",
    y = "Mean Ratio R",
    color = "Collinearity (rho)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_factor)

## 5) Stability plot: variance of R by method
var_plot_df <- replicate_results %>%
  mutate(
    rho = factor(rho, levels = sort(unique(rho))),
    pi = factor(pi, levels = sort(unique(pi))),
    beta_type = factor(beta_type, levels = c("sparse", "dense")),
    mechanism = factor(mechanism, levels = c("MCAR", "MAR"))
  ) %>%
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
  group_by(n, rho, pi, beta_type, mechanism, method) %>%
  summarise(
    R_var_mean = mean(R_var, na.rm = TRUE),
    .groups = "drop"
  )

p_var <- ggplot(
  var_plot_df,
  aes(x = pi, y = R_var_mean, color = method, group = method, shape = method)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  facet_grid(beta_type + rho ~ mechanism) +
  labs(
    title = "Stability of Relative MSE by Method",
    x = "Missingness proportion (pi)",
    y = "Mean Variance of R",
    color = "Method",
    shape = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

## 6) Observed vs Imputed Distribution (Mechanism Illustration)

library(mice)
library(dplyr)
library(ggplot2)

# Read representative missing-data example
dat <- readRDS("data/raw/missing_data_MAR_example.rds")

df <- as.data.frame(dat$X)
df$y <- dat$y

# Run one imputation for illustration
imp <- mice(df, m = 1, method = "pmm", printFlag = FALSE, seed = 123)
df_imp <- complete(imp)

# Choose one variable in the at-risk set
j <- 15

df_plot <- data.frame(
  value = c(df[[j]][!is.na(df[[j]])], df_imp[[j]][is.na(df[[j]])]),
  type = c(
    rep("Observed", sum(!is.na(df[[j]]))),
    rep("Imputed", sum(is.na(df[[j]])))
  )
)

mean_df <- df_plot %>%
  group_by(type) %>%
  summarise(mu = mean(value, na.rm = TRUE), .groups = "drop")

p_impdist <- ggplot(df_plot, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = mean_df,
    aes(xintercept = mu, color = type),
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = paste0("Observed vs Imputed Distribution for X", j),
    subtitle = "Illustration based on one representative imputed dataset",
    x = paste0("Value of X", j),
    y = "Density",
    fill = "Type",
    color = "Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_impdist)