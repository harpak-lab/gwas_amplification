# Here, I just want to test simulation of the new pgs corr methodology
set.seed(847)
'%>%' <- magrittr::'%>%'
library(ggplot2)

Sigma <- list()
#Sigma[[1]] <- (.25 ^ 2) * matrix(data = 1, nrow = 2, ncol = 2)
#Sigma[[2]] <- (.25 ^ 2) * matrix(data = c(1, 0, 0, 0), nrow = 2, ncol = 2)
#Sigma[[3]] <- (.25 ^ 2) * matrix(data = c(0, 0, 0, 1), ncol = 2)

#Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
Sigma[[1]] <- matrix(data = c(1, sqrt(1.5), sqrt(1.5), 1.5), ncol = 2)
#Sigma[[3]] <- matrix(data = c(1, sqrt(3), sqrt(3), 3), ncol = 2)

maf_simulator <- function(n) {

  pmax(pmin(rbeta(n, .5, 5), .5 - 1e-05), .01)

}

for (
  n in c(
    ceiling(1000 / .75),
    ceiling(2500 / .75),
    ceiling(5000 / .75),
    ceiling(10000 / .75),
    ceiling(15000 / .75),
    ceiling(25000 / .75),
    ceiling(50000 / .75),
    ceiling(100000 / .75)
  )
) {

  sim_out <- gamp::simulate_pgs_corr_one_env(
    n_e0 = n,
    n_e1 = n,
    n_snps = 500,
    maf_simulator,
    e0_h2 = .25,
    e1_h2 = .25,
    Sigma = Sigma,
    pi = c(1),
    num_sims = 15,
    s = .5,
    pval_thresh_additive = .1,
    pval_thresh_GxE = .1,
    selection_methods = c("additive", "GxE")
  )

  readr::write_rds(sim_out, glue::glue("rds_data/polymetric_{n}_zhu_test.rds"))

}

corr_vec <- c()
n_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()

for (selection_method in c("additive", "GxE")) {

  for (beta_method in c("additive", "GxE")) {

    for (n in c(
      ceiling(1000 / .75),
      ceiling(2500 / .75),
      ceiling(5000 / .75),
      ceiling(10000 / .75),
      ceiling(15000 / .75),
      ceiling(25000 / .75),
      ceiling(50000 / .75),
      ceiling(100000 / .75)
    )) {

      sim_data <- readr::read_rds(
        glue::glue("rds_data/polymetric_{n}_zhu_test.rds")
      )

      n_vec <- c(n_vec, n)
      selection_method_vec <- c(selection_method_vec, selection_method)
      beta_method_vec <- c(beta_method_vec, beta_method)
      corr_vec <- c(corr_vec, sim_data[[selection_method]][[beta_method]][['corr']])

    }

  }

}

corr_df <- data.frame(
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  n = n_vec,
  corr = corr_vec
)

corr_df <- corr_df %>%
  dplyr::filter(selection_method == beta_method)

ggplot(data =corr_df) +
  geom_point(aes(x = .75 * n, y = corr, color = selection_method), size = 1) +
  geom_line(aes(x = .75 * n, y = corr, color = selection_method)) +
  xlab("Number of Individuals in Training Set for Each Environment") +
  ylab("Environment 1 Test Set PGS Correlation") +
  labs(color="Model")
