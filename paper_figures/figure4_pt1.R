# code for figure 3 - analysis of pgs as % amplification increases
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)
library(foreach)

generate_fig3 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    # first, create list of covariance matrices
    Sigma <- list()
    Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
    Sigma[[2]] <- matrix(data = c(1, 0, 0, 0), nrow = 2, ncol = 2)
    Sigma[[3]] <- matrix(data = c(0, 0, 0, 1), ncol = 2)

    additive_cor_vec <- c()
    GxE_cor_vec <- c()
    mash_cor_vec <- c()
    additive_mse_vec <- c()
    GxE_mse_vec <- c()
    mash_mse_vec <- c()

    maf_simulator <- function(n) {

      pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), .005)

    }

    amp_range <- seq(0, 1, .05)
    len_amp_range <- length(amp_range)

    amp_corr_res <- foreach(
      i = 1:len_amp_range,
      .verbose = TRUE,
      .packages = c("gamp")
    ) %do% {

      pi <- c(1 - amp_range[i], amp_range[i] / 2, amp_range[i] / 2)

      gamp::simulate_pgs_corr_fast(
        n_e0 = 1.5e05,
        n_e1 = 1.5e05,
        n_snps = 1000,
        e0_h2 = .35,
        e1_h2 = .35,
        Sigma = Sigma,
        pi = pi,
        maf_simulator = maf_simulator,
        num_sims = 10,
        s = .5
      )

    }

    for (j in 1:length(amp_corr_res)) {

      additive_cor_vec <- c(additive_cor_vec, amp_corr_res[[j]]$additive_corr_pv4)
      GxE_cor_vec <- c(GxE_cor_vec, amp_corr_res[[j]]$GxE_corr_pv4)
      #mash_cor_vec <- c(mash_cor_vec, amp_corr_res[[j]]$mash_cor)
      additive_mse_vec <- c(additive_mse_vec, amp_corr_res[[j]]$additive_mse_pv4)
      GxE_mse_vec <- c(GxE_mse_vec, amp_corr_res[[j]]$GxE_mse_pv4)
      #mash_mse_vec <- c(mash_mse_vec, amp_corr_res[[j]]$mash_mse)

    }

    sim_df <- data.frame(
      amp = amp_range,
      additive_cor = additive_cor_vec,
      GxE_cor = GxE_cor_vec,
      #mash_cor = mash_cor_vec,
      additive_mse = additive_mse_vec,
      GxE_mse = GxE_mse_vec
      #mash_mse = mash_mse_vec
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig4_cond_spec.rds")

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
generate_fig3(write_rds = TRUE)

sim_df <- readr::read_rds("./rds_data/fig4_cond_spec.rds")

plot_df <- data.frame(
  pct_cond_spec = rep(sim_df$amp, 2),
  corr = c(sim_df$additive_cor, sim_df$GxE_cor),
  model_type = c(rep("additive", 21), rep("GxE", 21))
)


# res_df <- data.frame(
#   amp = c(rep(0, 2), rep(.1, 2), rep(.2, 2), rep(.3, 2), rep(.4, 2), rep(.5, 2)),
#   model_type = rep(c("Additive", "GxE"), 6),
#   corr = c(0.6268919, 0.616658, 0.6215744, 0.621064, 0.6150466, 0.6232702,
#            0.6131337, 0.6238186, 0.6099356, 0.6272873, 0.6051951, 0.6278665)
# )
#
# ggplot(data = sim_df, aes(x = amp, y = corr, color = model)) +
#  geom_point() +
#  geom_line()
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
# Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
# Sigma[[2]] <- matrix(data = c(1, sqrt(3), sqrt(3), 3), ncol = 2)
# Sigma[[3]] <- matrix(data = c(3, sqrt(3), sqrt(3), 1), ncol = 2)
# #Sigma[[4]] <- matrix(data = c(1, 0, 0, 0), ncol = 2)
# #Sigma[[5]] <- matrix(data = c(0, 0, 0, 1), ncol = 2)
# #Sigma[[6]] <- matrix(data = c(1, 0, 0, 1), ncol = 2)
# #Sigma[[7]] <- matrix(data = c(1, .5, .5, 1), ncol = 2)
#
#
# sim_out <- gamp::simulate_pgs_corr(
#   n_e0 = 10e04,
#   n_e1 = 10e04,
#   n_snps = 1000,
#   e0_h2 = .5,
#   e1_h2 = .5,
#   Sigma = Sigma,
#   pi = c((1/3), (1/3), (1/3)),
#   maf_simulator = maf_simulator,
#   num_sims = 10,
#   s = .5
# )

sim_df <- data.frame(
  amp = seq(0, 1, length.out = 7),
  additive_corr = c(0.5865399, 0.5854997, 0.5858993, 0.581313, 0.5801912, 0.5794931, 0.5783099),
  GxE_corr = c(0.5814596, 0.5826076, 0.5850098, 0.5821599, 0.5835542, 0.5827306, 0.5839103)
)

plot_df <- data.frame(
  amp = sim_df$amp,
  corr = c(sim_df$additive_corr, sim_df$GxE_corr),
  model_type = c(rep("additive", 7), rep("GxE", 7))
)

ggplot(data = plot_df, aes(x = pct_cond_spec, y = corr, color = model_type)) +
  geom_point() +
  geom_line()

