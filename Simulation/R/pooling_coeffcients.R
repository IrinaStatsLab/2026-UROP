# Simulation/R/pool_coeffcients.R
# Purpose: Pool coefficient vectors beta_1, ..., beta_m obtained from the imputed training datasets using different pooling strategies.

# Simple mean pooling 
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a pooled coefficient vector obtained by row-wise mean
pool_beta_mean <- function(beta_mat) {
  beta_pool <- rowMeans(beta_mat)
  return(beta_pool)
}

# SVD-based pooling
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a pooled coefficient vector based on the rank-1 SVD approximation of the coefficient matrix
pool_beta_svd <- function(beta_mat) {
  beta_mat <- as.matrix(beta_mat)
  sv <- svd(beta_mat)
  
  betamean <- rowMeans(beta_mat)
  v1 <- sv$v[, 1]
  
  beta_pool <- beta_mat %*% v1
  
  # sign alignment
  if (drop(crossprod(beta_pool, betamean)) < 0) {
    v1 <- -v1
    beta_pool <- -beta_pool
  }
  
  return(as.vector(beta_pool))
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