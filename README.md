# UROP 2026  
**Pooling Strategies for Ridge Regression with Multiple Imputation**

## Overview
This repository contains the implementation for a simulation study investigating how different pooling strategies affect the estimation of mean squared error (MSE) when ridge regression is combined with multiple imputation (MICE).

The primary objective is to evaluate how accurately different methods estimate MSE as a target parameter relative to the complete-data benchmark under missing data.

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

```text
R/
├── simulation_setting.R
├── data_generation.R
├── impose_missingness.R
├── analysis_utility.R
├── pooling_coefficients.R
├── one_simulation_replicate.R
├── train_test_replicate.R
├── split_data.R

analysis/
├── run_simulation.R
├── run_simulation_train_test.R
├── data_visualization.R

output/
├── figures/
├── rds/
```

---

## Simulation Frameworks

### 1. Baseline framework
- Script: `analysis/run_simulation.R`
- Purpose: Estimate MSE as a target parameter under imputed data
- Note: No train–test split is used

### 2. Train–test framework
- Script: `analysis/run_simulation_train_test.R`
- Purpose: Evaluate out-of-sample performance
- Note: Missingness is imposed on the training set only

---

## How to Run

Recommended order:

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

Simulation results are generated and saved during execution:

- Baseline results → `output/rds/replicate_results_*.rds`  
- Train–test results → `output/rds/replicate_results_train_test_*.rds`  

Plots are saved in:
- `output/figures/`

---

## Changing Simulation Settings

All settings are defined in:

```r
R/simulation_setting.R
```

Key parameters:

```r
n_grid       # sample size
pi_grid      # missingness proportion
rho_grid     # correlation level (collinearity)
beta_grid    # sparse vs dense
mech_grid    # MCAR vs MAR
m_imp        # number of imputations
R_reps       # number of replicates

lambda_mode  # "fixed_constant" or "complete_data_fixed"
lambda_fixed_value
```

---

## Important Notes

- **Baseline and train–test frameworks should be run separately**
- Changing any parameter (e.g., `pi_grid`, `rho_grid`, `mech_grid`, `lambda_mode`) requires rerunning the simulation
- Results are not automatically updated unless the simulation is rerun
- Parallel computing is used for efficiency