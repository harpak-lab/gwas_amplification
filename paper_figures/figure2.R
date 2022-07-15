# code for figure 2
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)

generate_fig2 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_diff_variance_2env(
      n_e0 = 1.5e05,
      n_e1 = 1.5e05,
      ref_se = .25,
      se_ratio_grid = seq(from = .5, to = 2, length.out = 100),
      fx_diff_norm_grid = seq(
        from = 0, to = 1, length.out = 100
      ) / gamp:::get_theo_obs_noise(se = .25, n = 1.5e05, maf = .4),
      sims_per_grid_pt = 1
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig2_df.rds")

  }

  sim_df <- sim_df %>%
    dplyr::mutate(mse_diff_e1 = mse_e1_additive - mse_e1_GxE)

  ggplot(data = sim_df, aes(x = fx_diff_norm, y = se_ratio, fill = mse_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    ylab(latex2exp::TeX("$\\frac{\\sigma_{1}}{\\sigma_{0}}$")) +
    xlab(latex2exp::TeX("$\\frac{|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|}{\\sigma_{0}}$"))

}

# first time simulation
generate_fig2(write_rds = TRUE)

# subsequent analysis
sim_df <- readr::read_rds("./rds_data/fig2_df.rds")
generate_fig2(sim_df)
