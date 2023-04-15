# Here, I'd like to create a simulation where the distribution of small effects is very
# different from the distribution of large effects

set.seed(8675309)
'%>%' <- magrittr::'%>%'

Sigma <- list()

Sigma[[1]] <- .02 * gamp:::make_amp_cov_mat(1, amp_coef = 3, amp_e0 = FALSE)
Sigma[[2]] <- .02 * gamp:::make_amp_cov_mat(1, amp_coef = 3, amp_e0 = TRUE)
#Sigma[[2]] <- .025 * matrix(data = c(1, -1, -1, 1), ncol = 2)
#Sigma[[3]] <- .025 * matrix(data = c(1, -.25, -.25, 1), ncol = 2)
#Sigma[[4]] <- .025 * matrix(data = c(0, 0, 0, 1), ncol = 2)
Sigma[[3]] <- matrix(data = c(1, 1, 1, 1), ncol = 2)

maf_simulator <- function(n) {
  
  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .001)
  
}

so <- gamp::simulate_perf_pval_par(
  n_e1 = 3333,
  n_e0 = 3333,
  n_snps = 10000,
  maf_simulator,
  e0_h2 = .2,
  e1_h2 = .2,
  Sigma = Sigma,
  pi = c(.95, .025, .025),
  num_sims = 10,
  s = .8,
  pval_thresh_vec = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-4, 1e-3, seq(.01, 1, .01))
)

n.cores <- 10
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

sim_out <- gamp::simulate_pgs_corr_fast_par(
  n_e0_train = 2500,
  n_e1_train = 2500,
  n_e0_test = 2500,
  n_e1_test = 2500,
  n_snps = 10000,
  maf_simulator,
  e0_h2 = .2,
  e1_h2 = .2,
  Sigma = Sigma,
  pi = c(.95, .025, .025),
  num_sims = 10,
  s = .8,
  pval_thresh_additive = .01,
  pval_thresh_GxE = .07
)

barplot(c(.332, .321, .406), names.arg = c("Additive", "GxE", "Mash"), ylim = c(0, .5))
