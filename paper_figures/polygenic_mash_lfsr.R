set.seed(8675309)
'%>%' <- magrittr::'%>%'
library(ggplot2)

Sigma <- list()
#Sigma[[1]] <- (.25 ^ 2) * matrix(data = 1, nrow = 2, ncol = 2)
#Sigma[[2]] <- (.25 ^ 2) * matrix(data = c(1, 0, 0, 0), nrow = 2, ncol = 2)
#Sigma[[3]] <- (.25 ^ 2) * matrix(data = c(0, 0, 0, 1), ncol = 2)

maf_simulator <- function(n) {

  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .005)

}

Sigma[[1]] <- (.25 ^ 2) * matrix(data = 1, nrow = 2, ncol = 2)
Sigma[[2]] <- (.25 ^ 2) * matrix(data = c(1, sqrt(1.25), sqrt(1.25), 1.25), ncol = 2)
Sigma[[3]] <- (.25 ^ 2) * matrix(data = c(1.25, sqrt(1.25), sqrt(1.25), 1), ncol = 2)


amp <- 1 / 3

sim_out <- gamp::simulate_mash_perf_lfsr(
  n_e0 = 13333,
  n_e1 = 13333,
  n_snps = 25000,
  maf_simulator,
  e0_h2 = .25,
  e1_h2 = .25,
  Sigma = Sigma,
  pi = c(1 - amp, amp/2, amp/2),
  num_sims = 25,
  s = .5
)

sim_out <- readr::read_rds("rds_data/mash_out_lfsr_preds_mse.rds")

thresh_vec <- c()
perc_selected_vec <- c()
corr_vec <- c()
mse_pred_vec <- c()

for (thresh in seq(.01, 1, .01)) {

  thresh_vec <- c(thresh_vec, thresh)
  perc_selected_vec <- c(perc_selected_vec, sim_out[[glue::glue("thresh_{thresh}")]][['perc_selected']])
  corr_vec <- c(corr_vec, sim_out[[glue::glue("thresh_{thresh}")]][['corr']])
  mse_pred_vec <- c(mse_pred_vec, sim_out[[glue::glue("thresh_{thresh}")]][['mse_preds']])

}

out_df <- data.frame(
  thresh = thresh_vec,
  perc_selected = perc_selected_vec,
  corr = corr_vec,
  mse_pred = mse_pred_vec
)

library(ggplot2)

ggplot(data = out_df, aes(x = thresh, y = mse_pred)) +
  geom_point() +
  geom_line() +
  xlab("LFSR Threshold") +
  ylab("Prediction MSE") +
  ggtitle("MASH")

readr::write_rds(sim_out, "rds_data/mash_out_lfsr_preds_mse.rds")
