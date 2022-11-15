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
    Sigma[[2]] <- matrix(data = c(1, sqrt(3), sqrt(3), 3), ncol = 2)

    additive_cor_vec <- c()
    GxE_cor_vec <- c()
    mash_cor_vec <- c()
    additive_mse_vec <- c()
    GxE_mse_vec <- c()
    mash_mse_vec <- c()

    maf_simulator <- function(n) {

      runif(n, .025, .5)

    }

    #amp_range <- seq(0, 1, .1)
    amp_range <- c(0, .1, .2, .3)
    len_amp_range <- length(amp_range)

    amp_corr_res <- foreach(
      i = 1:len_amp_range,
      .verbose = TRUE,
      .packages = c("gamp")
    ) %do% {

      pi <- c(1 - amp_range[i], amp_range[i])

      gamp::simulate_pgs_corr(
        n_e0 = 5000,
        n_e1 = 5000,
        n_snps = 1700,
        e0_h2 = .4,
        e1_h2 = .4,
        Sigma = Sigma,
        pi = pi,
        maf_simulator = maf_simulator,
        num_sims = 35,
      )

    }

    for (j in 1:length(amp_corr_res)) {

      additive_cor_vec <- c(additive_cor_vec, amp_corr_res[[j]]$additive_cor)
      GxE_cor_vec <- c(GxE_cor_vec, amp_corr_res[[j]]$GxE_cor)
      mash_cor_vec <- c(mash_cor_vec, amp_corr_res[[j]]$mash_cor)
      additive_mse_vec <- c(additive_mse_vec, amp_corr_res[[j]]$additive_mse)
      GxE_mse_vec <- c(GxE_mse_vec, amp_corr_res[[j]]$GxE_mse)
      mash_mse_vec <- c(mash_mse_vec, amp_corr_res[[j]]$mash_mse)

    }

    sim_df <- data.frame(
      amp = amp_range,
      additive_cor = additive_cor_vec,
      GxE_cor = GxE_cor_vec,
      mash_cor = mash_cor_vec,
      additive_mse = additive_mse_vec,
      GxE_mse = GxE_mse_vec,
      mash_mse = mash_mse_vec
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig3_df_pt1.rds")

  }

  # now, convert the dataframe to long format for ggplot
  # implement once I actually have the dataframe on hand


}

# first time simulation
generate_fig3(write_rds = TRUE)

# subsequent analysis
sim_df_1 <- readr::read_rds("./rds_data/fig3_df_pt1.rds")
sim_df_2 <- readr::read_rds("./rds_data/fig3_df_pt2.rds")
sim_df_3 <- readr::read_rds("./rds_data/fig3_df_pt3.rds")
sim_df <- rbind(sim_df_1, sim_df_2, sim_df_3)

generate_fig3(sim_df)


# temp code
plot_df <- data.frame(
  amp = rep(sim_df$amp, 3),
  corr = c(sim_df$additive_cor, sim_df$GxE_cor, sim_df$mash_cor),
  model_type = c(rep("additive", nrow(sim_df)), rep("GxE", nrow(sim_df)), rep("mash", nrow(sim_df)))
)

ggplot(data = plot_df, aes(x = amp, y = corr, color = model_type)) +
  geom_point() +
  geom_line()

# temp code
plot_df2 <- data.frame(
  amp = rep(sim_df$amp, 3),
  mse = c(sim_df$additive_mse, sim_df$GxE_mse, sim_df$mash_mse),
  model_type = c(rep("additive", nrow(sim_df)), rep("GxE", nrow(sim_df)), rep("mash", nrow(sim_df)))
)

ggplot(data = plot_df2, aes(x = amp, y = mse, color = model_type)) +
  geom_point() +
  geom_line()
