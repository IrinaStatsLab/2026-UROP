# R/save_example_data.R
# Purpose: generate and save representative datasets (complete + missing)

source("R/simulation_setting.R")
source("R/data_generation.R")
source("R/impose_missingness.R")

set.seed(123)

## Choose ONE scenario (representative)
n <- 1000
rho <- 0.5
pi <- 0.5
beta_type <- "sparse"
mechanism <- "MAR"

## Generate complete data
dat <- generate_complete_data(
  n = n, b = b, K = K, rho = rho,
  beta_type = beta_type,
  s = s, signal = signal,
  sigma_eps = sigma_eps,
  seed = 123
)

X_full <- dat$X
y <- dat$y

## Generate missing data
X_mis <- if (mechanism == "MCAR") {
  impose_mcar(X_full, pi, M_risk, seed = 456)
} else {
  impose_mar(X_full, y, pi, M_risk, seed = 789)
}

## Create folder
if (!dir.exists("data/raw")) {
  dir.create("data/raw", recursive = TRUE)
}

## Save
saveRDS(list(X = X_full, y = y),
        "data/raw/complete_data_example.rds")

saveRDS(list(X = X_mis, y = y),
        paste0("data/raw/missing_data_", mechanism, "_example.rds"))

cat("Example datasets saved in data/raw/\n")