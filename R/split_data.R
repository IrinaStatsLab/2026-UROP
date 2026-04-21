split_data <- function(complete_obj, train_prop = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.list(complete_obj) || !all(c("X", "y") %in% names(complete_obj))) {
    stop("Input must be output from generate_complete_data().")
  }
  
  dat <- as.data.frame(complete_obj$X)
  dat$y <- complete_obj$y
  
  n <- nrow(dat)
  idx <- sample(seq_len(n))
  
  n_train <- floor(train_prop * n)
  train_idx <- idx[1:n_train]
  test_idx  <- idx[(n_train + 1):n]
  
  list(
    train = dat[train_idx, , drop = FALSE],
    test  = dat[test_idx, , drop = FALSE],
    train_idx = train_idx,
    test_idx = test_idx
  )
}