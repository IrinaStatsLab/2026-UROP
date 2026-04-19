
# simulation/R/analysis_utility.R
ratios <- replicate_results$ratio_R

mu_hat  <- mean(ratios, na.rm = TRUE)
sd_hat  <- sd(ratios, na.rm = TRUE)

## Standardized ratios (CLT form)
Z_vals <- (ratios - mu_hat) / sd_hat


hist(Z_vals, breaks = 40, freq = FALSE,
     main = "Standardized Ratio R",
     xlab = "Z")

curve(dnorm(x, 0, 1), add = TRUE, col = "red", lwd = 2)


#variance MC
var_imp <- var(mse_imp, na.rm = TRUE)

delta_var_ratio <- var_imp / (length(na.omit(mse_imp)) * mse_full^2)

delta_sd_ratio <- sqrt(delta_var_ratio)


library(dplyr)
library(ggplot2)
library(readr)

df1 <- MAR_and_DENSE %>%
  mutate(mechanism = "MAR", beta_type = "Dense")

df2 <- MAR_and_Sparse %>%
  mutate(mechanism = "MAR", beta_type = "Sparse")

df3 <- MCAR_and_DENSE %>%
  mutate(mechanism = "MCAR", beta_type = "Dense")

df4 <- MCAR_and_Sparse %>%
  mutate(mechanism = "MCAR", beta_type = "Sparse")

df <- bind_rows(df1, df2, df3, df4)

ggplot(df, aes(x = factor(pi), y = ratio_R_mean)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_grid(beta_type ~ mechanism) +
  labs(
    x = "Missingness proportion (π)",
    y = "R = MSE_pooled / MSE_full"
  ) +
  theme_minimal(base_size = 14)

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

df1 <- read_csv("results_dense_MAR.csv") %>%
  mutate(mechanism = "MAR", beta_type = "Dense")

df2 <- read_csv("results_sparse_MAR.csv") %>%
  mutate(mechanism = "MAR", beta_type = "Sparse")

df3 <- read_csv("results_table_dense_MCAR.csv") %>%
  mutate(mechanism = "MCAR", beta_type = "Dense")

df4 <- read_csv("results_sparse_MCAR.csv") %>%
  mutate(mechanism = "MCAR", beta_type = "Sparse")

df <- bind_rows(df1, df2, df3, df4) %>%
  mutate(pi = as.numeric(pi))

df_long <- df %>%
  select(pi, mechanism, beta_type, ratio_R_mean, ratio_R_coef_mean) %>%
  pivot_longer(
    cols = c(ratio_R_mean, ratio_R_coef_mean),
    names_to = "method",
    values_to = "R"
  ) %>%
  mutate(
    method = recode(
      method,
      ratio_R_mean = "Pool MSE",
      ratio_R_coef_mean = "Pool Coefficient"
    )
  )

ggplot(df_long, aes(x = pi, y = R, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(beta_type ~ mechanism) +
  labs(
    x = "Missingness proportion (π)",
    y = "R",
    color = "Pooling Strategy"
  ) +
  theme_minimal(base_size = 14)
theme(
  legend.position = "bottom",
)
  