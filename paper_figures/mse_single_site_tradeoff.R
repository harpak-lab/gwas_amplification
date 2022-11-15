'%>%' <- magrittr::'%>%'
library(ggplot2)

generate_mse_tradeoff_single_site_plot <- function(
  sim_df = NULL, write_rds = FALSE, sigma_ratio = 1
) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_var_ratio_2env(
      sigma_ratio = sigma_ratio,
      se_grid = seq(
        from = .001, to = 1, length.out = 100
      ),
      fx_diff_grid = seq(from = 0, to = 1, by = .025)
    )

  }

  if(write_rds) {

    readr::write_rds(
      sim_df, glue::glue("./rds_data/mse_tradeoff_single_site_{sigma_ratio}.rds")
    )

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff_e0 = mse_e0_additive - mse_e0_GxE,
      bias_diff_e0 = bias_e0_additive ^ 2 - bias_e0_GxE ^ 2,
      var_diff_e0 = var_e0_additive - var_e0_GxE
    )

  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Difference in Model Mean Squared Error") +
    geom_segment(
      aes(x = .015, y = 0.01, xend = 1, yend = .71), linetype = "dashed"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="MSE Diff.") +
    annotate("text", x = .35, y = .75, label = "Additive Better", size = 5) +
    annotate("text", x = .75, y = .15, label = "GxE Better", size = 5)


  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = bias_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle(latex2exp::TeX("Difference in Model Bias")) +
    labs(fill="Bias Diff.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate(
      "text", x = .225, y = .5,
      label = latex2exp::TeX("$Bias(Add) \\approx Bias(GxE)$"), size = 2.3
      ) +
    annotate("text", x = .775, y = .5, label = "Bias(Add) > Bias(GxE)", size = 2.4)

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = var_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Difference in Model Variance") +
    scale_fill_gradient2() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Var Diff.") +
    annotate("text", x = .5, y = .9, label = "Var(GxE) > Var(Add)", size = 5) +
    annotate("text", x = .5, y = .1,
             label = latex2exp::TeX("$Var(Add) \\approx Var(GxE)$"), size = 5)


  return(
    list(
      mse_plot = g1,
      bias_plot = g2,
      variance_plot = g3
    )
  )

}

# first time simulation
plot <- generate_mse_tradeoff_single_site_plot(write_rds = T, sigma_ratio = 1)


# subsequent analysis
sigma_ratio <- 1
sim_df <- readr::read_rds("./rds_data/mse_tradeoff_single_site_{sigma_ratio}.rds")

g <- generate_fig1(sim_df)
