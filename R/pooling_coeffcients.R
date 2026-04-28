
# R/pool_coeffcients.R
# Pooling strategies for combining coefficient estimates across imputations

# Mean pooling of coefficients
pool_beta_mean <- function(beta_mat) {
  beta_pool <- rowMeans(beta_mat)
  return(beta_pool)
}

# SVD-based pooling of coefficients (captures dominant direction across imputations)
beta_pool_svd <- function(beta_mat) {
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
  
  # normalized weights (preserve relative importance)
  w <- v1 / sum(abs(v1))
  
  # weighted pooling
  beta_pool <- beta_mat %*% w
  
  # align sign with mean direction
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

# Wrapper for pooling methods
pool_beta_all <- function(beta_mat) {
  return(list(
    mean = pool_beta_mean(beta_mat),
    svd  = pool_beta_svd(beta_mat)
  ))
}