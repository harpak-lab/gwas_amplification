# Here, I'd like to create a simulation where the distribution of small effects is very
# different from the distribution of large effects

set.seed(8675309)
'%>%' <- magrittr::'%>%'

Sigma <- list()

Sigma[[1]] <- matrix(data = c(1, sqrt(1.5), sqrt(1.5), 1.5), ncol = 2)
#Sigma[[2]] <- .025 * matrix(data = c(1, -1, -1, 1), ncol = 2)
#Sigma[[3]] <- .025 * matrix(data = c(1, -.25, -.25, 1), ncol = 2)
#Sigma[[4]] <- .025 * matrix(data = c(0, 0, 0, 1), ncol = 2)
Sigma[[2]] <- .1 * matrix(data = c(1, 1, 1, 1), ncol = 2)

maf_simulator <- function(n) {
  
  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .001)
  
}

sim_out <- gamp::simulate_pgs_corr_fast(
  n_e0_train = 2500,
  n_e1_train = 2500,
  n_e0_test = 2500,
  n_e1_test = 2500,
  n_snps = 10000,
  maf_simulator,
  e0_h2 = (.4 ^ 2),
  e1_h2 = (.35 ^ 2),
  Sigma = Sigma,
  pi = c(.01, .99),
  num_sims = 10,
  s = .8,
  pval_thresh_additive = .1,
  pval_thresh_GxE = .1
)
