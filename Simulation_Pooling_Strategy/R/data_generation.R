# Simulation_Pooling_Strategy/R/data_generation.R
# Purpose: block auto regression covariance construction, true coefficient generation, and data generation.

## Logistic inverse function (Used in MAR missingness later)
inv_logit <- function(x) 1 / (1 + exp(-x))

## Generate block auto regression covariance matrix
## Within each block, each predictor has auto regression correlation with the others; 
## Across blocks, each predictor is independent of the others.
## b   : block size
## K   : number of blocks
## rho : within-block auto regression correlation parameter (0 < rho < 1)
## Returns: p x p covariance matrix, where p = b * K

make_sigma_block_ar1 <- function(b, K, rho) {
  Sigma_block <- outer(1:b, 1:b, function(i, j) rho^abs(i - j))
  Sigma <- matrix(0, nrow = b*K, ncol = b*K)
  for (k in 1:K) {
    idx <- ((k-1)*b + 1):(k*b)
    Sigma[idx, idx] <- Sigma_block
  }
  Sigma
}

## Construct true coefficient vector beta
##
## p         : number of predictors
## beta_type : "sparse" or "dense"
## s         : number of non-zero coefficients when sparse
## signal    : signal magnitude
##
## Returns: beta (length p) true coefficient vector

make_beta <- function(p, beta_type = c("sparse", "dense"),
                      s = 10, signal = 1.0, seed = NULL) {
  beta_type <- match.arg(beta_type)
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  beta <- rep(0, p)
  if (beta_type == "sparse") {
    idx <- 1:s
    signs <- sample(c(-1, 1), size = s, replace = TRUE)
    beta[idx] <- signs * signal
  } else {
    beta[] <- signal / sqrt(p)
  }
  beta
}

## Generate complete data under linear model y = X beta + eps
## setting must contain: n, b, K, rho, beta_type, s, signal, sigma_eps

generate_complete_data <- function(n, b, K, rho, beta_type,
                                   s = 10, signal = 1.0,
                                   sigma_eps = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  p <- b * K
  Sigma <- make_sigma_block_ar1(b, K, rho)
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", seq_len(p))
  
  beta <- make_beta(p, beta_type, s, signal,
                    seed = if (is.null(seed)) NULL else seed + 999)
  
  y <- as.numeric(X %*% beta + rnorm(n, mean = 0, sd = sigma_eps))
  data <- data.frame(y = y, X)
  list(X = X, y = y, beta = beta, Sigma = Sigma, data = data)
}