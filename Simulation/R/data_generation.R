
# simulation/R/data_generation.R
# Contains: AR covariance, beta construction, and data generation.

round_int <- function(x) as.integer(round(x))

# Generate AR covariance matrix for predictors
# p: number of predictors
# rho: correlation strength (0 < rho < 1)
make_sigma_ar1 <- function(p, rho) {
  outer(1:p, 1:p, function(i, j) rho^abs(i - j))
}

# Construct true coefficient vector beta
# beta_type: "sparse" or "dense"
# s: number of non-zero coefficients when sparse
# signal: magnitude of signal
make_beta <- function(p,
                      beta_type = c("sparse", "dense"),
                      s = 10,
                      signal = 1.0,
                      seed = NULL) {
  beta_type <- match.arg(beta_type)
  
  if (!is.null(seed)) set.seed(seed)
  
  beta <- rep(0, p)
  
  if (beta_type == "sparse") {
    idx <- seq_len(min(s, p))
    beta[idx] <- signal * sample(c(-1, 1), length(idx), replace = TRUE)
  } else {
    beta <- rnorm(p, mean = 0, sd = signal / sqrt(p))
  }
  
  beta
}

# Generate complete data (X, y) under linear model y = X beta + eps
# setting must contain: n, p, rho, beta_type, s, signal, sigma_eps
generate_data <- function(setting, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- setting$n
  p <- setting$p
  rho <- setting$rho
  
  Sigma <- make_sigma_ar1(p, rho)
  
  # MASS::mvrnorm generates multivariate normal predictors
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  
  beta <- make_beta(
    p = p,
    beta_type = setting$beta_type,
    s = setting$s,
    signal = setting$signal,
    seed = if (is.null(seed)) NULL else seed + 100
  )
  
  y <- as.numeric(X %*% beta + rnorm(n, sd = setting$sigma_eps))
  
  list(X_full = X, y = y, beta_true = beta)
}

