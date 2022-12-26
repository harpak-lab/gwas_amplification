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

  color_scale_limits <- c(
    min(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0)),
    max(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0))
  )
  
  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("MSE of Environment 0 Estimate") +
    geom_segment(
      aes(x = .015, y = 0.01, xend = 1, yend = .71), linetype = "dashed"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Difference (au)") +
    annotate("text", x = .35, y = .75, label = "Additive Better", size = 5) +
    annotate("text", x = .75, y = .15, label = "GxE Better", size = 5) +
    scale_fill_gradient2(limits = color_scale_limits)
  
  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = bias_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle(latex2exp::TeX("Bias of Environment 0 Estimate")) +
    labs(fill="Bias Diff.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = .165, y = .54, label = "equal", size = 4) +
    annotate("text", x = .165, y = .46, label = "bias", size = 4) +
    annotate("text", x = .77, y = .54, label = "higher bias", size = 4) +
    annotate("text", x = .77, y = .46, label = "with additive", size = 4) +
    geom_segment(aes(x = .6, y = .35, xend = .95, yend = .35),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    scale_fill_gradient2(limits = color_scale_limits)

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = var_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Variance of Environment 0 Estimate") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Var Diff.") +
    annotate("text", x = .5, y = .9, label = "higher variance", size = 4) +
    annotate("text", x = .5, y = .82, label = "with GxE", size = 4) +
    annotate("text", x = .5, y = .1,
             label = "equal variance", size = 4) +
    geom_segment(aes(x = .13, y = .75, xend = .13, yend = .95),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    scale_fill_gradient2(limits = color_scale_limits)
  
  return(
    list(
      mse_plot = g1,
      bias_plot = g2,
      variance_plot = g3
    )
  )

}

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
  
  color_scale_limits <- c(
    min(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0)),
    max(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0))
  )
  
  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = color_scale_limits) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("MSE of Environment 0 Estimate") +
    geom_curve(
      aes(x = 0, y = 1.74, xend = 1.5, yend = .78),
      curvature = -.24,
      linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Difference (au)") +
    annotate("text", x = .5, y = 1, label = "Additive Better", size = 5) +
    annotate("text", x = 1, y = 2, label = "GxE Better", size = 5) + 
    theme(aspect.ratio=1)
  
  
  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = bias_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = color_scale_limits) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Bias of Enviroment 0 Estimate") +
    labs(fill="Bias Diff.") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(aspect.ratio=1) +
    annotate("text", x = .165, y = 1.35, label = "equal", size = 4) +
    annotate("text", x = .165, y = 1.2, label = "bias", size = 4) +
    annotate("text", x = 1.2, y = 1.35, label = "higher bias", size = 4) +
    annotate("text", x = 1.2, y = 1.2, label = "with additive", size = 4) +
    geom_segment(aes(x = .97, y = 1, xend = 1.45, yend = 1),
                 arrow = arrow(length = unit(0.5, "cm")))
  
  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = var_diff_e0)) +
    geom_tile() +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Variance of Environment 0 Estimate") +
    scale_fill_gradient2(limits = color_scale_limits) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Var Diff.") +
    annotate("text", x = .75, y = 2.25, label = "lower variance", size = 4) +
    annotate("text", x = .75, y = 2.05, label = "with GxE", size = 4) +
    annotate("text", x = .75, y = .85,
             label = "higher variance", size = 4) +
    annotate("text", x = .75, y = .65,
             label = "with GxE", size = 4) +
    geom_segment(aes(x = .3, y = .95, xend = .3, yend = .55),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    geom_segment(aes(x = .3, y = 1.95, xend = .3, yend = 2.35),
                 arrow = arrow(length = unit(0.5, "cm"))) +
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
plot1 <- generate_mse_tradeoff_single_site_se_ratio_plot(write_rds = T)

# first time simulation
plot2 <- generate_mse_tradeoff_single_site_plot(write_rds = T, sigma_ratio = 1)

library(ggpubr)
library(grid)

fig <- ggarrange(
  plot2$mse_plot + rremove("ylab") + rremove("xlab"),
  plot2$bias_plot + rremove("ylab") + rremove("xlab"), 
  plot2$variance_plot + rremove("ylab") + rremove("xlab"),
  nrow = 1,
  common.legend = TRUE,
  legend = "right",
  labels = "AUTO"
)

annotate_figure(fig, 
                left = textGrob("Environment Specific Standard Deviation (au)", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 9)),
                bottom = textGrob("Difference in Environment Specific Effects (au)", gp = gpar(cex = 1, fontsize = 9), hjust = .62))


fig2 <- ggarrange(
  plot1$mse_plot + rremove("ylab") + rremove("xlab"),
  plot1$bias_plot + rremove("ylab") + rremove("xlab"), 
  plot1$variance_plot + rremove("ylab") + rremove("xlab"),
  nrow = 1,
  common.legend = TRUE,
  legend = "right",
  labels = "AUTO"
)

annotate_figure(fig2, 
                left = textGrob("(Environment 0 SD) / (Environment 1 SD)", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 9)),
                bottom = textGrob("Difference in Environment Specific Effects (au)", gp = gpar(cex = 1, fontsize = 9), hjust = .62))
