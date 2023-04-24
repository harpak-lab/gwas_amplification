# Here, I want to make a lower power and high power version of the figure
# I think that I can do 2x amplification
# but I don't want to hard code anything
# I can build a script that takes in the following arguments
# (1) number of individuals
# (2) number of snps
# (3) heritability
# (4) quantity of amplification
# (5) percentage of amplification
# (6) Number of simulations

command_args = commandArgs(trailingOnly = TRUE)
n_individuals = as.integer(command_args[1])
n_snps = as.integer(command_args[2])
heritability = as.double(command_args[3])
amp_coef = as.double(command_args[4])
amp_pct = as.double(command_args[5])
n_sims = as.integer(command_args[6])

set.seed(1)
'%>%' <- magrittr::'%>%'

Sigma <- list()
Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
Sigma[[2]] <- (1 / amp_coef) * gamp:::make_amp_cov_mat(
  desired_corr = 1, 
  amp_coef = amp_coef
)

maf_simulator <- function(n) {
  
  .5 * rbeta(n, 1, 5)
  
}

sim_out_list <- gamp::simulate_pgs_corr_fast(
    n_e0_train = n_individuals,
    n_e1_train = n_individuals,
    n_e0_test = 3000,
    n_e1_test = 3000,
    n_snps = n_snps,
    maf_simulator,
    e0_h2 = heritability,
    e1_h2 = heritability,
    Sigma = Sigma,
    pi = c(1 - amp_pct, amp_pct),
    num_sims = n_sims,
    n_ascertained_snps = floor(n_snps * .5 * .5),
    s = .5,
    selection_methods = c("additive", "GxE"),
    beta_methods = c("additive", "GxE")
)
  
readr::write_rds(
  sim_out_list, 
  glue::glue("pgs_sim_n_{n_individuals}_n_snps_{n_snps}_h2",
             "_{heritability}_amp_coef_{amp_coef}_amp_perc_",
             "{amp_pct}_n_sims_{n_sims}.rds")
)
