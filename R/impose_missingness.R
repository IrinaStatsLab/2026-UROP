
# R/impose_missingness.R
# Missing-data mechanisms used in the simulation (MCAR and MAR)

# Impose MCAR missingness: entries are missing independently with probability pi
impose_mcar <- function(X, pi, M_risk, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xmis <- X
  for (j in M_risk) {
    miss_idx <- runif(nrow(Xmis)) < pi
    Xmis[miss_idx, j] <- NA
  }
  Xmis
}

# Impose MAR missingness: missingness depends on observed covariates and outcome
impose_mar <- function(X, y, pi, M_risk,
                       a1 = 0.5, a2 = 0.5, a3 = 0.5,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xmis <- X
  
  # Fully observed covariates used to drive MAR mechanism
  x1 <- X[, 1]
  x2 <- X[, 2]
  
  # Calibrate intercept to achieve target missingness level pi
  f_root <- function(alpha0) {
    pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
    mean(1 - pr_obs) - pi
  }
  
  alpha0 <- tryCatch(
    uniroot(f_root, interval = c(-20, 20))$root,
    error = function(e) NA_real_
  )
  
  if (is.na(alpha0)) alpha0 <- 0
  
  # Observation probability under MAR
  pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
  
  # Apply missingness to selected predictors
  for (j in M_risk) {
    miss_idx <- runif(length(y)) > pr_obs
    Xmis[miss_idx, j] <- NA
  }
  Xmis
}