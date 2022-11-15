# code for figure 2
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)

generate_fig8 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_diff_maf_2env(
      n_e0 = 1.5e03,
      n_e1 = 1.5e03,
      maf_grid = seq(from = .01, to = .1, by = .0001),
      fx_diff_grid = seq(from = 0, to = .5, by = .01),
      sims_per_grid_pt = 5,
      sigma_e0 = 7.5,
      sigma_e1 = 15
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig8_df.rds")

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff_e0 = mse_e0_additive - mse_e0_GxE,
      var_diff_e0 = var_e0_additive - var_e0_GxE,
      bias_diff_e0 = bias_e0_additive ^ 2 - bias_e0_GxE,
      mse_diff_e1 = mse_e1_additive - mse_e1_GxE,
      var_diff_e1 = var_e1_additive - var_e1_GxE,
      bias_diff_e1 = bias_e1_additive ^ 2 - bias_e1_GxE,
    )

  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = maf, fill = mse_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    ggtitle("MSE Additive - MSE GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  g2 <- ggplot(data = sim_df, aes(x = fx_diff_norm, y = se_ratio, fill = bias_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    ylab(latex2exp::TeX("$\\frac{\\sigma_{1}}{\\sigma_{0}}$")) +
    xlab(latex2exp::TeX("$\\frac{|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|}{\\sigma_{0}}$")) +
    ggtitle("Bias^2 Additive - Bias^2 GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  g3 <- ggplot(data = sim_df, aes(x = fx_diff_norm, y = se_ratio, fill = var_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    ylab(latex2exp::TeX("$\\frac{\\sigma_{1}}{\\sigma_{0}}$")) +
    xlab(latex2exp::TeX("$\\frac{|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|}{\\sigma_{0}}$")) +
    ggtitle("Var Additive - Var GxE") +
    theme(plot.title = element_text(hjust = 0.5))


}

# first time simulation
generate_fig8(write_rds = TRUE)

# subsequent analysis
sim_df <- readr::read_rds("./rds_data/fig8_df.rds")
generate_fig8(sim_df)
