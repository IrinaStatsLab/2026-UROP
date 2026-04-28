
# R/simulation_setting.R
# Global settings for simulation study

# Load required packages
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

# Data-generating setup
b <- 10                         # block size
K <- 5                          # number of blocks
p <- b * K                      # total number of predictors (50)
M_risk <- 11:30                 # 20 predictors can be missing; others are anchors

# Scenario grids
beta_grid <- c("dense")         # "sparse" / "dense"
mech_grid <- c("MAR")           # "MAR" / "MCAR"
rho_grid <- c(0.2, 0.5, 0.8)    # collinearity levels
pi_grid  <- c(0.2, 0.5, 0.8)    # missingness proportion
n_grid   <- c(500, 1000)        # sample size

# Model Settings
s <- 10                         # number of nonzero coefficients in sparse setting
signal <- 1.0                   # signal strength
sigma_eps <- 1.0                # noise SD in outcome model
lambda_mode <- "fixed_constant" #"complete_data_fixed" / "fixed_constant"
lambda_fixed_value <- 1.0

# Multiple imputation settings
m_imp <- 5                    
R_reps <- 200                 

# Simulation control
seed0 <- 123
ncores <- max(1, parallel::detectCores() - 1)


