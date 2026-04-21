# R/split_data.R
# Purpose: Split the complete dataset into training and test sets.

split_data <- function(complete_obj, train_prop = 0.8, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.list(complete_obj) || !all(c("X", "y") %in% names(complete_obj))) {
    stop("Input must be output from generate_complete_data().")
  }
  
  dat <- as.data.frame(complete_obj$X)
  colnames(dat) <- colnames(complete_obj$X)
  dat$y <- complete_obj$y
  
  split <- initial_split(dat, prop = train_prop)
  
  list(
    train = training(split),
    test  = testing(split)
  )
}