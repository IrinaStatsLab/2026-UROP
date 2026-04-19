
# R/impose_missingness.R
# Impose missingness on selected predictors
#
# MCAR (Missing Completely At Random): Each entry is missing independently of all variables
#
# Settings:
# X       : n x p complete predictor matrix
# pi      : target missing proportion (0 < pi < 1)
# M_risk  : indices of predictors subject to missingness
# seed    : optional seed for reproducibility
#
# Returns: Xmis, which denotes the predictor matrix with MCAR missingness imposed

impose_mcar <- function(X, pi, M_risk, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xmis <- X
  # For each at-risk predictor, independently mask entries with prob = pi
  for (j in M_risk) {
    miss_idx <- runif(nrow(Xmis)) < pi
    Xmis[miss_idx, j] <- NA
  }
  Xmis
}

# MAR (Missing At Random): Missingness depends on observed variables (x1, x2, y), but not on the unobserved value itself.
#
# Settings:
# X       : n x p complete predictor matrix
# y       : response vector
# pi      : target missing proportion (0 < pi < 1)
# M_risk  : indices of predictors subject to missingness
# a1,a2,a3: coefficients controlling MAR dependence
#
# Returns: Xmis, which denotes the predictor matrix with MAR missingness imposed

impose_mar <- function(X, y, pi, M_risk,
                       a1 = 0.5, a2 = 0.5, a3 = 0.5,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xmis <- X
  
  # Fully observed covariates used to drive MAR mechanism
  x1 <- X[, 1]
  x2 <- X[, 2]
  
  # Solve once for intercept so that E[missing] ≈ pi
  f_root <- function(alpha0) {
    pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
    mean(1 - pr_obs) - pi
  }
  
  alpha0 <- tryCatch(
    uniroot(f_root, interval = c(-20, 20))$root,
    error = function(e) NA_real_
  )
  
  if (is.na(alpha0)) alpha0 <- 0
  
  # Shared observation probability under MAR
  pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
  
  # For each at-risk predictor, independently generate missingness
  for (j in M_risk) {
    miss_idx <- runif(length(y)) > pr_obs
    Xmis[miss_idx, j] <- NA
  }
  
  Xmis
}