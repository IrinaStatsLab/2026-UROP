
# R/simulation_setting.R

# Install Missing Packages and Load All Required Packages
## Notes on key packages:
## - MASS: multivariate normal generator (mvrnorm) used in the data generating mechanism
## - mice: multiple imputation (mice)
## - glmnet: ridge fitting (glmnet) and optional CV tuning (cv.glmnet)
## - dplyr: result wrangling (bind_rows, group_by, summarize, etc.)
## - foreach + doParallel: parallel simulation replications
## - doRNG:guarantees reproducible random number streams when running simulations in parallel.
## - rsample: train test split
pkgs <- c("MASS", "mice", "glmnet", "dplyr", "foreach", 
          "doParallel", "doRNG", "rsample")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(MASS) 
library(mice)
library(glmnet) 
library(dplyr) 
library(foreach) 
library(doParallel) 
library(doRNG) 
library(rsample)

## Global constants (Protocol)
b <- 10                         # block size
K <- 5                          # number of blocks
p <- b * K                      # total number of predictors (50)
M_risk <- 11:30                 # 20 predictors can be missing; others are anchors

# Fixed Scenario grids
beta_grid <- c("sparse")
mech_grid <- c("MAR")

# Varied Scenario grids
rho_grid <- c(0.2, 0.5, 0.8) # collinearity levels
pi_grid  <- c(0.2, 0.5, 0.8) # missingness proportion among M_risk variables
n_grid   <- c(500, 1000) # sample size

## Settings
s <- 10                         # number of nonzero coefficients in sparse setting
signal <- 1.0                   # coefficient magnitude / signal strength
sigma_eps <- 1.0                # noise SD in outcome model
lambda_mode <- "complete_data_fixed"
lambda_fixed_value <- 1.0

## Feature selection settings
perm_m <- 30                     # permutations per replicate (null max|beta|)
sel_alpha <- 0.95                # quantile for selection threshold

## MI settings
m_imp <- 5                    
R_reps <- 200                   # number of Monte Carlo replications 

## Simulation settings
seed0 <- 123
ncores <- max(1, parallel::detectCores() - 1)


