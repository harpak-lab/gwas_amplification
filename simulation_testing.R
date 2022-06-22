source("simulation.R")

Sigma <- list()
Sigma[[1]] <- diag(2)
Sigma[[2]] <- matrix(data = c(1, 1, 1, 2), nrow = 2, ncol = 2, byrow = TRUE)

sim_out <- simulate_2env_gwas_mse(
  n_sims = 2,
  n_snps = 1000,
  n_e0 = 1000,
  n_e1 = 1000,
  sigma_e0 = 1,
  sigma_e1 = 1,
  Sigma = Sigma,
  pi = NULL,
  s = 0.8,
  corr_range = seq(from = -1, to = 1, by = .5),
  amp_range = c(1.5, 2),
  int_e0 = 0,
  int_e1 = 0,
  maf_e0 = .4,
  maf_e1 = .4
)
