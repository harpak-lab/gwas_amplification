n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

library(foreach)
library(dplyr)

maf_simulator <- function(n) {
  
  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .01)
  
}

Sigma_init <- list()
Sigma_init[[1]] <- matrix(data = c(1, 1, 1, 1), ncol = 2)

no_amp_sim <- gamp::simulate_pgs_corr_fast_par(
  n_e0_train = 2500,
  n_e1_train = 2500,
  n_e0_test = 2500,
  n_e1_test = 2500,
  n_snps = 10000,
  maf_simulator,
  e0_h2 = .4,
  e1_h2 = .4,
  Sigma = Sigma_init,
  pi = c(1),
  num_sims = 128,
  s = .8,
  pval_thresh_additive = 5e-3,
  pval_thresh_GxE = .07,
  beta_methods = c("additive", "GxE"),
  selection_methods = c("additive")
)

perc_amp_vec <- c(seq(0, 1, length.out = 50), rep(0, 50))
amp_size_vec <- c(rep(1, 50), seq(1, 3, length.out = 50))
additive_corr_vec <- c(rep(no_amp_sim$corr[1], 100))
GxE_corr_vec <- c(rep(no_amp_sim$corr[2], 100))

for (perc_amp in seq(0, 1, length.out = 50)[-1]) {
  
  print(glue::glue("Working on {perc_amp} amplification."))
  
  for (amp in seq(1, 3, length.out = 50)[-1]) {
    
    Sigma <- list()
    Sigma[[1]] <- matrix(data = c(1, 1, 1, 1), ncol = 2)
    Sigma[[2]] <- matrix(data = c(1, sqrt(amp), sqrt(amp), amp), ncol = 2)
    pi <- c(1 - perc_amp, perc_amp)
    
    sim_out <- gamp::simulate_pgs_corr_fast_par(
      n_e0_train = 2500,
      n_e1_train = 2500,
      n_e0_test = 2500,
      n_e1_test = 2500,
      n_snps = 10000,
      maf_simulator,
      e0_h2 = .4,
      e1_h2 = .4,
      Sigma = Sigma,
      pi = pi,
      num_sims = 128,
      s = .8,
      pval_thresh_additive = 5e-3,
      pval_thresh_GxE = .07,
      beta_methods = c("additive", "GxE"),
      selection_methods = c("additive")
    )
    
    perc_amp_vec <- c(perc_amp_vec, perc_amp)
    amp_size_vec <- c(amp_size_vec, amp)
    additive_corr_vec <- c(additive_corr_vec, sim_out$corr[1])
    GxE_corr_vec <- c(GxE_corr_vec, sim_out$corr[2])
    
  }
  
}

sim_res_df <- data.frame(
  perc_amp = perc_amp_vec,
  amp = amp_size_vec,
  additive_corr = additive_corr_vec,
  GxE_corr = GxE_corr_vec
)

sim_res_df <- sim_res_df %>%
  mutate(additive_improvement = (additive_corr - GxE_corr) / GxE_corr)

readr::write_csv(sim_res_df, "./colorscale_add_vs_GxE.csv")
