#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

# Here, I want to do some sort of an ascertainment study

# code for figure 3 - analysis of pgs as % amplification increases
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)
library(foreach)

generate_fig6 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    # first, create list of covariance matrices

    Sigma <- list()

    Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
    Sigma[[2]] <- matrix(data = c(1, 0, 0, 0), nrow = 2, ncol = 2)
    Sigma[[3]] <- matrix(data = c(0, 0, 0, 1), ncol = 2)

    additive_power_vec <- c()
    GxE_power_vec <- c()

    additive_tp_mse_vec <- c()
    GxE_tp_mse_vec <- c()

    maf_simulator <- function(n) {

      pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .005)

    }

    amp_range <- seq(from = 0, to = 1, by = .05)
    len_amp_range <- length(amp_range)

    thresh_asc_res <- foreach(
      i = 1:len_amp_range,
      .verbose = TRUE,
      .packages = c("gamp")
    ) %do% {

      gamp::simulate_ascertainment_fast(
        n_e0 = 1.5e05,
        n_e1 = 1.5e05,
        n_snps = 1000,
        e0_h2 = .35,
        e1_h2 = .35,
        Sigma = Sigma,
        pi = c(1 - amp_range[i], amp_range[i] / 2, amp_range[i] / 2),
        maf_simulator = maf_simulator,
        num_sims = 100,
        pval_thresh = 5e-8,
        s = .5
      )

    }

    for (j in 1:length(thresh_asc_res)) {

      additive_power_vec <- c(additive_power_vec, thresh_asc_res[[j]]$additive_power)
      GxE_power_vec <- c(GxE_power_vec, thresh_asc_res[[j]]$GxE_power)

      additive_tp_mse_vec <- c(additive_tp_mse_vec, thresh_asc_res[[j]]$additive_tp_mse)
      GxE_tp_mse_vec <- c(GxE_tp_mse_vec, thresh_asc_res[[j]]$GxE_tp_mse)


    }

    sim_df <- data.frame(
      amp = amp_range,

      additive_power = additive_power_vec,
      additive_tp_mse = additive_tp_mse_vec,

      GxE_power = GxE_power_vec,
      GxE_tp_mse = GxE_tp_mse_vec
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, glue::glue("./rds_data/fig6_ind.rds"))

  }

  # now, convert the dataframe to long format for ggplot
  # implement once I actually have the dataframe on hand


}

#df1 <- readr::read_rds("./rds_data/fig4_pt1.rds")
#df2 <- readr::read_rds("./rds_data/fig4_pt2.rds")

#df <- rbind(df1, df2)

#sim_df <- data.frame(
#  amp = rep(sim$amp, 3),
#  model = c(rep("Additive", 11), rep("GxE", 11), rep("Mash", 11)),
#  corr = c(sim$additive_cor, sim$GxE_cor, sim$mash_cor)
#)

# first time simulation
generate_fig6(write_rds = TRUE)

sim_df <- readr::read_rds("./rds_data/fig6_ind.rds")

power_plot_df <- data.frame(
  pct_cond_spec = rep(sim_df$amp, 2),
  power = c(sim_df$additive_power, sim_df$GxE_power),
  model_type = c(rep("additive", 21), rep("GxE", 21))
)

ggplot(power_plot_df) +
  geom_point(aes(x = pct_cond_spec, y = power, color = model_type)) +
  geom_line(aes(x = pct_cond_spec, y = power, color = model_type))

mse_plot_df <- data.frame(
  pct_cond_spec = rep(sim_df$amp, 2),
  mse = c(sim_df$additive_tp_mse, sim_df$GxE_tp_mse),
  model_type = c(rep("additive", 21), rep("GxE", 21))
)

ggplot(mse_plot_df) +
  geom_point(aes(x = pct_cond_spec, y = mse, color = model_type)) +
  geom_line(aes(x = pct_cond_spec, y = mse, color = model_type))

# sim_df <- readr::read_rds("./rds_data/fig6.rds")
#
# plot_df <- data.frame(
#   pval = rep(sim_df$pval_thresh, 2),
#   power = c(sim_df$additive_asc_tpr, sim_df$GxE_asc_tpr),
#   model = c(rep("additive", 4), rep("GxE", 4))
# )
#
# ggplot(data = plot_df) +
#   geom_point(aes(x = -log10(pval), y = power, color = model)) +
#   geom_line(aes(x = -log10(pval), y = power, color = model))
#
#
# maf_simulator <- function(n) {
#
#   pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .005)
#
# }
#
# Sigma <- list()
#
# #Sigma[[1]] <- .05 * matrix(data = 1, nrow = 2, ncol = 2)
# Sigma[[1]] <- .05 * matrix(data = c(1, sqrt(2), sqrt(2), 2), ncol = 2)
#
#
# pgs_out <- gamp::simulate_pgs_corr(
#   n_e0 = 1.5e04,
#   n_e1 = 1.5e04,
#   n_snps = 1700,
#   e0_h2 = .6,
#   e1_h2 = .6,
#   Sigma = Sigma,
#   pi = c(1),
#   maf_simulator = maf_simulator,
#   num_sims = 10,
#   s = .75
# )







sim_df <- data.frame(
  amp = seq(0, 1, length.out = 7),
  additive_power = c(0.5598, 0.5578, 0.5564, 0.554, 0.5574, 0.5522, 0.5502),
  GxE_power = c(0.4412, 0.4356, 0.4369, 0.4333, 0.4321, 0.4309, 0.4274)
)

sim_df2 <- data.frame(
  amp = seq(0, 1, length.out = 7),
  additive_mse = c(0.01194095, 0.02379697, 0.03887729, 0.05097351, 0.06271446, 0.07677285, 0.08899077),
  GxE_mse = c(0.02062767, 0.02063835, 0.02356353, 0.02395724, 0.02579837, 0.0269602, 0.0297422)
)

plot_df <- data.frame(
  amp = sim_df$amp,
  power = c(sim_df$additive_power, sim_df$GxE_power),
  model_type = c(rep("additive", 7), rep("GxE", 7))
)

ggplot(data = plot_df, aes(x = amp, y = power, color = model_type)) +
  geom_point() +
  geom_line()

plot_df2 <- data.frame(
  amp = sim_df2$amp,
  mse = c(sim_df2$additive_mse, sim_df2$GxE_mse),
  model_type = c(rep("additive", 7), rep("GxE", 7))
)

ggplot(data = plot_df2, aes(x = amp, y = mse, color = model_type)) +
  geom_point() +
  geom_line()
