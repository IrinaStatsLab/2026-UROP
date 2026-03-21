#  Simulation_Pooling_Strategy/R/run_mice_imputation.R
## Purpose: Run multiple imputation by chained equations (MICE) on the incomplete training data and return the completed datasets for ridge regression fitting.

# Run MICE on incomplete training data
# data_mis : incomplete training data frame with y and predictors
# m        : number of imputations
# maxit    : number of MICE iterations
# method   : optional imputation method specification
# seed     : optional seed for reproducibility
# Returns: A list containing the original mids object returned by mice() and a list of completed datasets for each imputation.

run_mice_imputation <- function(data_mis, m, maxit, method = NULL, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # If no custom method is supplied, let mice choose defaults
  if (is.null(method)) {
    imp <- mice::mice(
      data = data_mis,
      m = m,
      maxit = maxit,
      printFlag = FALSE,
      seed = seed
    )
  } else {
    imp <- mice::mice(
      data = data_mis,
      m = m,
      maxit = maxit,
      method = method,
      printFlag = FALSE,
      seed = seed
    )
  }
  
  # Extract each completed dataset
  imputed_list <- lapply(seq_len(m), function(j) {
    mice::complete(imp, action = j)
  })
  
  return(list(
    mids = imp,
    imputed_list = imputed_list
  ))
}