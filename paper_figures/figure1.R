# code for figure 1
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)

generate_fig1 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_fx_diff_noise_tradeoff_2env(
      n_e0 = 1.5e05,
      n_e1 = 1.5e05,
      se_grid = seq(from = .025, to = 1, by = .025),
      fx_diff_grid = seq(from = 0, to = 1, by = .025),
      sims_per_grid_pt = 1
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig1_df.rds")

  }

  sim_df <- sim_df %>%
    dplyr::mutate(mse_diff_e1 = mse_e1_additive - mse_e1_GxE)

  ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    geom_segment(aes(x = .030265, y = 0.02500095, xend = 1, yend = 0.683035), linetype = "dashed") +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$"))

}

# first time simulation
generate_fig1(write_rds = TRUE)

# subsequent analysis
sim_df <- readr::read_rds("./rds_data/fig1_df.rds")
generate_fig1(sim_df)
