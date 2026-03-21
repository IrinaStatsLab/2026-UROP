# Simulation_Pooling_Strategy/R/split_data.R
# Purpose: Split the complete dataset into training, validation, and test sets.

## Split complete data into training / validation / test sets
## complete_obj : the full list returned by generate_complete_data(), or the data frame containing y and predictors
## train_prop   : proportion of observations assigned to training set
## val_prop     : proportion assigned to validation set
## seed         : optional random seed
## Returns: A list containing train(training data frame), validation(validation data frame), and test(test data frame).

split_data <- function(complete_obj, train_prop, val_prop, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract data frame
  if (!is.list(complete_obj) || !("data" %in% names(complete_obj))) {
    stop("Input must be output from generate_complete_data().")
  }
  dat <- complete_obj$data
  
  # Random split the complete dataset into three disjoint subsets
  # To ensure the split is unbiased and independent of ordering
  n <- nrow(dat)
  shuffled_idx <- sample(seq_len(n))
  
  # Number of observations in each subset
  n_train <- floor(train_prop * n)
  n_val   <- floor(val_prop * n)
  
  # Assign indices to each subset
  train_idx <- shuffled_idx[1:n_train]
  validation_idx <- shuffled_idx[(n_train + 1):(n_train + n_val)]
  test_idx <- shuffled_idx[(n_train + n_val + 1):n]
  
  # Extract corresponding data subsets
  train_data <- dat[train_idx, , drop = FALSE]
  validation_data <- dat[validation_idx, , drop = FALSE]
  test_data <- dat[test_idx, , drop = FALSE]
  
  # Return split results
  return(list(
    train = train_data,
    validation = validation_data,
    test = test_data,
    train_idx = train_idx,
    validation_idx = validation_idx,
    test_idx = test_idx
  ))
}