
# Simulation_Pooling_Strategy/R/impose_missingness.R
# Impose missingness on selected predictors

# MCAR (Missing Completely At Random): Each entry is missing independently of all variables
# Settings:
# data   : complete training data frame with y and predictors
# pi     : target missing proportion (0 < pi < 1)
# M_risk : indices of predictors subject to missingness
# seed   : optional seed for reproducibility
# Returns: A data frame with MCAR missingness imposed on selected predictors

impose_mcar <- function(data, pi, M_risk, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat_mis <- data
  
  # Predictor column names corresponding to M_risk
  x_names <- paste0("X", M_risk)
  
  if (!all(x_names %in% colnames(dat_mis))) {
    stop("Some predictors in M_risk are not found in the data.")
  }
  
  # For each selected predictor, independently impose missingness with probability pi for each observation
  for (xj in x_names) {
    miss_idx <- runif(nrow(dat_mis)) < pi
    dat_mis[miss_idx, xj] <- NA
  }
  return(dat_mis)
}

# MAR (Missing At Random): Missingness depends on observed variables (x1, x2, y), but not on the unobserved value itself.
#
# Settings:
# data   : complete training data frame with y and predictors
# pi     : target missing proportion (0 < pi < 1)
# M_risk : indices of predictors subject to missingness
# a1,a2,a3: coefficients controlling MAR dependence of missingness on X1, X2, and y, respectively
#
# Returns: A data frame with MAR missingness imposed on selected predictors

impose_mar <- function(data, pi, M_risk,
                       a1 = 0.5, a2 = 0.5, a3 = 0.5,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat_mis <- data
  
  x_names <- paste0("X", M_risk)
  
  x1 <- dat_mis$X1
  x2 <- dat_mis$X2
  y  <- dat_mis$y
  
  for (xj in x_names) {
    # Solve for intercept alpha0 so that the expected missing proportion is approximately pi.
    # pr_obs denotes the probability that an entry remains observed.
    f_root <- function(alpha0) {
      pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
      mean(1 - pr_obs) - pi
    }
    
    alpha0 <- tryCatch(uniroot(f_root, interval = c(-20, 20))$root,
                       error = function(e) NA_real_)
    # Fallback if root-finding fails 
    if (is.na(alpha0)) alpha0 <- 0
    # Observation probability under MAR
    pr_obs <- inv_logit(alpha0 + a1 * x1 + a2 * x2 + a3 * y)
    # Generate missing indicators
    miss_idx <- runif(nrow(dat_mis)) > pr_obs
    dat_mis[miss_idx, xj] <- NA
  }
  return(dat_mis)
}

# Wrapper function for missingness mechanism "MCAR" or "MAR"
impose_missingness <- function(data, mechanism, pi, M_risk,
                               a1 = NULL, a2 = NULL, a3 = NULL,
                               seed = NULL) {
  if (mechanism == "MCAR") {
    return(impose_mcar(data = data,pi = pi,M_risk = M_risk,seed = seed))
  }
  
  if (mechanism == "MAR") {
    if (is.null(a1) || is.null(a2) || is.null(a3)) {
      stop("For MAR, a1, a2, and a3 must be provided.")
    }
    
    return(impose_mar(data = data, pi = pi, M_risk = M_risk, a1 = a1, a2 = a2, a3 = a3, seed = seed))
  }
}