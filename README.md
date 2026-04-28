# UROP 2026  
**Pooling Strategies for Ridge Regression with Multiple Imputation**

## Overview

This repository contains the implementation for a simulation study investigating how different pooling strategies affect the estimation of mean squared error (MSE) when ridge regression is combined with multiple imputation (MICE).

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
├── tables/
```

---

## Simulation Frameworks

### 1. Baseline framework
- Script: `analysis/run_simulation.R`
- Purpose: Estimate MSE as a target parameter under imputed data  
- Pipeline:
  - Uses `data_generation.R`, `impose_missingness.R`, `analysis_utility.R`, and `pooling_coefficients.R`
  - Uses `one_simulation_replicate.R` for each simulation replicate  
  - Does **not** use `split_data.R` or `train_test_replicate.R`  
- Note:
  - All required functions are automatically sourced within the script

---

### 2. Train–test framework
- Script: `analysis/run_simulation_train_test.R`
- Purpose: Evaluate out-of-sample performance  
- Pipeline:
  - Uses `data_generation.R`, `impose_missingness.R`, `analysis_utility.R`, and `pooling_coefficients.R`
  - Uses `split_data.R` and `train_test_replicate.R`  
  - Does **not** use `one_simulation_replicate.R`  
- Note:
  - All required functions are automatically sourced within the script

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
n_grid
pi_grid
rho_grid
beta_grid
mech_grid
m_imp
R_reps

lambda_mode
lambda_fixed_value
```

---

## Important Notes

- Baseline and train–test frameworks should be run separately  
- Changing any simulation setting requires rerunning the simulation  
- Results will not update automatically unless the simulation is rerun  
- Parallel computing is used for efficiency  
