# Simulation_Pooling_Strategy/R/pool_beta.R
# Purpose: Pool coefficient vectors beta_1, ..., beta_m obtained from the imputed training datasets using different pooling strategies.

# Simple mean pooling 
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a pooled coefficient vector obtained by row-wise mean
pool_beta_mean <- function(beta_mat) {
  beta_pool <- rowMeans(beta_mat)
  return(beta_pool)
}

# Geometric mean pooling
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# eps      : small positive constant to avoid log(0)
# Returns a pooled coefficient vector obtained by applying a sign-preserving geometric mean row-wise
pool_beta_geo <- function(beta_mat, eps = 1e-8) {
  beta_pool <- apply(beta_mat, 1, function(b) {
    s <- sign(mean(b))
    mag <- exp(mean(log(abs(b) + eps)))
    s * mag
  })
  return(beta_pool)
}

# SVD-based pooling
# beta_mat : coefficient matrix with rows = coefficients and columns = imputations
# Returns a pooled coefficient vector based on the rank-1 SVD approximation of the coefficient matrix
pool_beta_svd <- function(beta_mat) {
  sv <- svd(beta_mat)
  
  # coefficient space
  u1 <- sv$u[, 1]
  # variation size
  d1 <- sv$d[1]
  # weights for every imputation
  v1 <- sv$v[, 1]
  
  # Rank-1 reconstruction
  beta_rank1 <- d1 * u1 %*% t(v1)
  
  # Average reconstructed coefficients across imputations
  beta_pool <- rowMeans(beta_rank1)
  names(beta_pool) <- rownames(beta_mat)
  
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
    geo  = pool_beta_geo(beta_mat),
    svd  = pool_beta_svd(beta_mat)
  ))
}