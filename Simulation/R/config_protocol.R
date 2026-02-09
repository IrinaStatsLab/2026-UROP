
# Install Missing Packages and Load All Required Packages
## Notes on key packages:
## - MASS: multivariate normal generator (mvrnorm) used in the data generating mechanism
## - mice: multiple imputation (mice), MAR/MCAR missingness via ampute, completed data via complete
## - glmnet: ridge fitting (glmnet) and optional CV tuning (cv.glmnet)
## - dplyr: result wrangling (bind_rows, group_by, summarize, etc.)
## - foreach + doParallel: parallel simulation replications
## - doRNG:guarantees reproducible random number streams when running simulations in parallel.
pkgs <- c("MASS", "mice", "glmnet", "dplyr", "foreach", "doParallel", "doRNG")
for (pkg in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(pkg)
}
library(MASS) 
library(mice)
library(glmnet) 
library(dplyr) 
library(foreach) 
library(doParallel) 
library(doRNG) 

## Global constants (Protocol)
n <- 1000                       # sample size for the main (analysis) dataset
b <- 10                         # block size
K <- 5                          # number of blocks
p <- b * K                      # total number of predictors (50)
M_risk <- 11:30                 # 20 predictors can be missing; others are anchors

## Scenario grids (factorial design)
rho_grid <- c(0.2, 0.5, 0.8)    # Autoregression correlation strength among predictors
pi_grid  <- c(0.2, 0.5, 0.8)    # missingness proportion among M_risk variables
beta_grid <- c("sparse", "dense")
mech_grid <- c("MCAR", "MAR")

## Settings
s <- 10                         # number of nonzero coefficients in sparse setting
signal <- 1.0                   # coefficient magnitude / signal strength
sigma_eps <- 1.0                # noise SD in outcome model
lambda_fixed <- 1.0             # Ridge penalty 

## MI settings
m_imp <- 20                     # number of imputations (use 5-10 for debugging)
R_reps <- 200                   # number of Monte Carlo replications (use 1000 for final runs)

## Prediction Evaluation Setting:
n_test <- 1000                  # sample size for the independent test set

