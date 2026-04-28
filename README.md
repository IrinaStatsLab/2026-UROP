# UROP 2026  
**Pooling Strategies for Ridge Regression with Multiple Imputation**

## Overview
This repository contains the implementation for a simulation study investigating how different pooling strategies affect the estimation of mean squared error (MSE) when ridge regression is combined with multiple imputation (MICE).

The focus is on estimating MSE as a target parameter under missing data.

Two evaluation frameworks are implemented:
- **Baseline framework**: evaluates how well different methods estimate MSE on imputed datasets relative to the complete-data benchmark
- **Train–test framework**: evaluates out-of-sample performance using a data split

---

## Dependencies
All required R packages will be automatically installed if missing.

Main packages:
- MASS
- mice
- glmnet
- dplyr
- foreach
- doParallel
- doRNG
- rsample

---

## Code Organization

The code is organized as follows:

```text
R/
├── simulation_setting.R        # global simulation parameters
├── data_generation.R           # generate X, beta, y
├── impose_missingness.R        # MCAR / MAR missingness mechanisms
├── analysis_utility.R          # ridge regression, MSE, diagnostics, MICE wrapper
├── pooling_coefficients.R      # pooling methods (mean, SVD)
├── one_simulation_replicate.R  # one baseline simulation replicate
├── train_test_replicate.R      # one train-test simulation replicate
├── split_data.R                # train-test split function

analysis/
├── run_simulation.R            # run baseline simulation
├── run_simulation_train_test.R # run train-test simulation
├── data_visualization.R        # generate figures and tables

output/
├── figures/                    # plots
├── rds/                        # replicate-level results
├── tables/                     # summary tables
```

---

## How to Run

Run the scripts in the following order.

### Step 1: Baseline simulation
```r
source("analysis/run_simulation.R")
```

### Step 2: Train–test simulation
```r
source("analysis/run_simulation_train_test.R")
```

### Step 3: Visualization
```r
source("analysis/data_visualization.R")
```

---

## Results

All results are saved in the `output/` folder:

- `output/rds/`: replicate-level simulation results  
- `output/tables/`: summarized results across replicates  
- `output/figures/`: plots used for analysis  

Raw result files can be inspected using:

```r
readRDS("output/rds/filename.rds")
```

---

## Changing Simulation Settings

All simulation settings are defined in:

```r
R/simulation_setting.R
```

Key parameters you may want to modify:

```r
# sample size
n_grid <- c(500, 1000)

# missingness level
pi_grid <- c(0.2, 0.5, 0.8)

# correlation
rho_grid <- c(0.2, 0.5, 0.8)

# coefficient structure
beta_grid <- c("sparse", "dense")

# missingness mechanism
mech_grid <- c("MCAR", "MAR")

# number of imputations
m_imp <- 5

# number of simulation replicates
R_reps <- 200

# lambda setting
lambda_mode <- "fixed_constant"   # or "complete_data_fixed"
lambda_fixed_value <- 1.0
```

---

## Notes
- Parallel computing is used for efficiency.
