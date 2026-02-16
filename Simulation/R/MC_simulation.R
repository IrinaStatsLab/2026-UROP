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
