'%>%' <- magrittr::'%>%'
library(ggplot2)

generate_mse_tradeoff_single_site_se_ratio_plot <- function(
  sim_df = NULL, write_rds = FALSE
) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_var_ratio_grid_2env(
      sigma_ratio_grid = seq(
        from = .4, to = 2.5, length.out = 150
      ),
      fx_diff_grid = seq(from = 0, to = 1.5, by = .025)
    )

  }

  if(write_rds) {

    readr::write_rds(
      sim_df, glue::glue("./rds_data/mse_tradeoff_single_site_sig_ratio.rds")
    )

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff_e0 = mse_e0_additive - mse_e0_GxE,
      bias_diff_e0 = bias_e0_additive ^ 2 - bias_e0_GxE ^ 2,
      var_diff_e0 = var_e0_additive - var_e0_GxE,
      mse_diff_e1 = mse_e1_additive - mse_e1_GxE,
      bias_diff_e1 = bias_e1_additive ^ 2 - bias_e1_GxE ^ 2,
      var_diff_e1 = var_e1_additive - var_e1_GxE
    )

  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Difference in Model Mean Squared Error") +
    geom_curve(
      aes(x = 0, y = 1.74, xend = 1.5, yend = .78),
      curvature = -.24,
      linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="MSE Diff.") +
    annotate("text", x = .5, y = 1, label = "Additive Better", size = 5) +
    annotate("text", x = 1, y = 2, label = "GxE Better", size = 5) + 
    theme(aspect.ratio=1)


  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = bias_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Difference in Model Bias") +
    labs(fill="Bias Diff.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(
      "text", x = .25, y = 1.5,
      label = latex2exp::TeX("$Bias(Add) \\approx Bias(GxE)$"), size = 2.3
    ) +
    annotate("text", x = 1.25, y = 1.5, label = "Bias(Add) > Bias(GxE)", size = 2.4) +
    theme(aspect.ratio=1)

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = var_diff_e0)) +
    geom_tile() +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Difference in Model Variance") +
    scale_fill_gradient2() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Var Diff.") +
    annotate("text", x = .75, y = 2.25, label = "Var(GxE) > Var(Add)", size = 5) +
    annotate("text", x = .75, y = .75,
             label = "Var(GxE) < Var(Add)", size = 5) +
    theme(aspect.ratio=1)


  return(
    list(
      mse_plot = g1,
      bias_plot = g2,
      variance_plot = g3
    )
  )

}

# first time simulation
plot <- generate_mse_tradeoff_single_site_se_ratio_plot(write_rds = T)

ggarrange(plot$mse_plot, plot$bias_plot, plot$variance_plot)


# subsequent analysis
sigma_ratio <- 1
sim_df <- readr::read_rds("./rds_data/mse_tradeoff_single_site_{sigma_ratio}.rds")

g <- generate_fig1(sim_df)
