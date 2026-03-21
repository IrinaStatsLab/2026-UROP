# Simulation_Pooling_Strategy/R/make_setting.R
# Purpose: Create a setting list for one scenario and one replicate

make_setting <- function(n, rho, pi, beta_type, mechanism, rep_id = 1) {
  list(n = n, b = b, K = K, rho = rho, beta_type = beta_type, s = s, signal = signal, sigma_eps = sigma_eps,
       train_prop = train_prop, val_prop = val_prop, mechanism = mechanism, pi = pi, M_risk = M_risk,
       a1 = a1, a2 = a2, a3 = a3, m = m, maxit = maxit, imputation_method = imputation_method,
       lambda_grid = lambda_grid,
       
    # Seed To Ensure Each Replicate is Different
    seed_data = seed0 + 10000 * rep_id + 1,
    seed_split = seed0 + 10000 * rep_id + 2,
    seed_missing = seed0 + 10000 * rep_id + 3,
    seed_mice = seed0 + 10000 * rep_id + 4)
}