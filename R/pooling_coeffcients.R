# R/pool_coeffcients.R
# Purpose: Pool coefficient vectors beta_1, ..., beta_m obtained from the imputed training datasets using different pooling strategies.

# Simple mean pooling 
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a pooled coefficient vector obtained by row-wise mean
pool_beta_mean <- function(beta_mat) {
  beta_pool <- rowMeans(beta_mat)
  return(beta_pool)
}

# SVD-based pooling of regression coefficients across imputations
# beta_mat: p x m matrix
#   - rows    = coefficients/features
#   - columns = imputations
# returns:
#   - pooled coefficient vector of length p
pool_beta_svd_weighted_signed <- function(beta_mat) {
  beta_mat <- as.matrix(beta_mat)
  
  if (!is.numeric(beta_mat)) {
    stop("beta_mat must be a numeric matrix.")
  }
  
  if (ncol(beta_mat) < 2) {
    stop("beta_mat must contain at least two imputations (columns).")
  }
  
  # SVD: B = UDV^T
  sv <- svd(beta_mat)
  
  # first right singular vector
  v1 <- sv$v[, 1]
  
  # mean coefficient vector for sign alignment
  beta_bar <- rowMeans(beta_mat)
  
  # signed normalized weights (preserve sign, control scale)
  w <- v1 / sum(abs(v1))
  
  # weighted pooling
  beta_pool <- beta_mat %*% w
  
  # optional sign alignment for interpretability
  if (drop(crossprod(beta_pool, beta_bar)) < 0) {
    beta_pool <- -beta_pool
    w <- -w
  }
  
  beta_pool <- as.vector(beta_pool)
  
  if (!is.null(rownames(beta_mat))) {
    names(beta_pool) <- rownames(beta_mat)
  }
  
  return(beta_pool)
}

# Wrapper function to pool coefficients using all methods
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a list containing pooled coefficient vectors from:
# mean : simple mean pooling
# geo  : geometric mean pooling
# svd  : SVD-based pooling
pool_beta_all <- function(beta_mat) {
  return(list(
    mean = pool_beta_mean(beta_mat),
    svd  = pool_beta_svd(beta_mat)
  ))
}