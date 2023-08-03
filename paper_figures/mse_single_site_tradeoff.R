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

  #color_scale_limits <- c(
  #  min(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0)),
  #  max(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0))
  #)
  
  color_scale_limits <- c(
    min(min(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0)), -.5),
    max(c(sim_df$mse_diff_e0, sim_df$bias_diff_e0, sim_df$var_diff_e0))
  )
  
  slope <- 1 / (2 * sqrt(.75 - (1 / (4 * sigma_ratio))))
  
  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Mean Squared Error (MSE)") +
    geom_segment(
      aes(x = .015, y = 0.015, xend = 1, yend = 1 * slope), linetype = "dashed"
    )  +
    labs(fill="Difference") +
    annotate("text", x = .4, y = .7, label = "additive better", size = 7) +
    annotate("text", x = .75, y = .15, label = "GxE better", size = 7) +
    scale_fill_gradient2(limits = color_scale_limits) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14), 
      legend.title=element_text(size=14),
      legend.title.align=0.5
    )
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = bias_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle(latex2exp::TeX("Bias")) +
    labs(fill="Difference") +
    annotate("text", x = .165, y = .55, label = "equal", size = 7) +
    annotate("text", x = .165, y = .45, label = "bias", size = 7) +
    annotate("text", x = .725, y = .55, label = "higher bias", size = 7) +
    annotate("text", x = .725, y = .45, label = "with additive", size = 7) +
    geom_segment(aes(x = .575, y = .325, xend = .925, yend = .325),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    scale_fill_gradient2(limits = color_scale_limits) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14), 
      legend.title=element_text(size=14),
      legend.title.align=0.5
    )
  

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = var_diff_e0)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific Standard Error") +
    xlab("Difference in Environment Specific Effects") +
    ggtitle("Variance") +
    labs(fill="Difference") +
    annotate("text", x = .5, y = .75, label = "higher variance", size = 7) +
    annotate("text", x = .5, y = .65, label = "with GxE", size = 7) +
    annotate("text", x = .5, y = .1,
             label = "equal variance", size = 5.5) +
    geom_segment(aes(x = .09, y = .55, xend = .09, yend = .8),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    scale_fill_gradient2(limits = color_scale_limits) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14), 
      legend.title=element_text(size = 14),
      legend.title.align=0.5
    )
  
  return(
    list(
      mse_plot = g1,
      bias_plot = g2,
      variance_plot = g3
    )
  )

}

generate_multi_ratio_plot <- function() {
  
  sim_df_r1 <- gamp:::sim_var_ratio_2env_v2(
    r = 1,
    se_grid = seq(
      from = .001, to = 1, length.out = 100
    ),
    fx_diff_grid = seq(from = 0, to = 1, by = .025)
  )
  
  sim_df_r1 <- sim_df_r1 %>%
    dplyr::mutate(
      mse_diff_e0 = mse_additive - mse_GxE
    )
  
  sim_df_r_two_thirds <- gamp:::sim_var_ratio_2env_v2(
    r = .5,
    se_grid = seq(
      from = .001, to = 1, length.out = 100
    ),
    fx_diff_grid = seq(from = 0, to = 1, by = .025)
  )
  
  sim_df_r_two_thirds <- sim_df_r_two_thirds %>%
    dplyr::mutate(
      mse_diff_e0 = mse_additive - mse_GxE
    )
  
  sim_df_r_one_pt_five <- gamp:::sim_var_ratio_2env_v2(
    r = 2,
    se_grid = seq(
      from = .001, to = 1, length.out = 100
    ),
    fx_diff_grid = seq(from = 0, to = 1, by = .025)
  )
  
  sim_df_r_one_pt_five <- sim_df_r_one_pt_five %>%
    dplyr::mutate(
      mse_diff_e0 = mse_additive - mse_GxE
    )
  
  #color_scale_limits <- c(
  #  min(c(sim_df_r1$mse_diff_e0, sim_df_r_half$mse_diff_e0, sim_df_r2$mse_diff_e0)),
  #  max(c(sim_df_r1$mse_diff_e0, sim_df_r_half$mse_diff_e0, sim_df_r2$mse_diff_e0))
  #)
  
  g1 <- ggplot(data = sim_df_r1, aes(x = fx_diff, y = se_A, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-.65, .25),
                         breaks = seq(-.65, .25, .15),
                         labels = seq(-.65, .25, .15)) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Equal Estimation Noise") +
    annotate("text", x = .4, y = .7, label = "additive better", size = 5.5) +
    annotate("text", x = .75, y = .15, label = "GxE better", size = 5.5) +
    geom_segment(
      aes(x = .001, y = 0.001, xend = 1, yend = 1/sqrt(2)), linetype = "dashed"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(fill="MSE Difference\n(additive - GxE)") +
    theme(aspect.ratio=1, legend.title=element_text(size=14))
  
  g2 <- ggplot(data = sim_df_r_one_pt_five, aes(x = fx_diff, y = se_A, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-.65, .25),
                         breaks = seq(-.65, .25, .15),
                         labels = seq(-.65, .25, .15)) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Higher Estimation Noise in Focal Context") +
    geom_segment(
      aes(x = .001, y = 0.001, xend = 1, yend = 1 / sqrt(3 - (1/2))), linetype = "dashed"
    ) +
    annotate("text", x = .425, y = .65, label = "additive better", size = 5.5) +
    annotate("text", x = .75, y = .15, label = "GxE better", size = 5.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(fill="MSE Difference\n(additive - GxE)") +
    theme(aspect.ratio=1, legend.title=element_text(size=14))
  
  g3 <- ggplot(data = sim_df_r_two_thirds, aes(x = fx_diff, y = se_A, fill = mse_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-.65, .25),
                         breaks = seq(-.65, .25, .15),
                         labels = seq(-.65, .25, .15)) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Lower Estimation Noise in Focal Context") +
    geom_segment(
      aes(x = .001, y = 0.001, xend = 1, yend = 1 / sqrt(3 - (1/.5))), linetype = "dashed"
    ) +
    annotate("text", x = .35, y = .7, label = "additive better", size = 5.5) +
    annotate("text", x = .75, y = .15, label = "GxE better", size = 5.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(fill="MSE Difference\n(additive - GxE)") +
    theme(aspect.ratio=1, legend.title=element_text(size=14))
  
  return(
    list(
      mse_plot_r1 = g1,
      mse_plot_rone_five = g2,
      mse_plot_rtwo_thirds = g3
    )
  )
  
}

generate_mse_tradeoff_single_site_se_ratio_plot <- function(
    sim_df = NULL, write_rds = FALSE
) {
  
  if(is.null(sim_df)) {
    
    sim_df <- gamp::sim_var_ratio_grid_2env(
      sigma_ratio_grid = seq(
        from = .75, to = 2.25, length.out = 150
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
    scale_fill_gradient2(limits = c(-.75, 1.125),
                         breaks = seq(-.75, 1.125, .375),
                         labels = seq(-.75, 1.125, .375)) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("MSE") +
    geom_curve(
      aes(x = 0, y = 1.74, xend = 1.5, yend = .78),
      curvature = -.3,
      linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Difference") +
    annotate("text", x = .5, y = 1, label = "additive better", size = 6) +
    annotate("text", x = 1.05, y = 1.85, label = "GxE better", size = 6) + 
    theme(aspect.ratio=1, legend.title=element_text(size=8.5))
  
  
  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = bias_diff_e0)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-.75, 1.125),
                         breaks = seq(-.75, 1.125, .375),
                         labels = seq(-.75, 1.125, .375)) +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Bias") +
    labs(fill="Difference") +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = .165, y = 1.56, label = "equal", size = 6) +
    annotate("text", x = .165, y = 1.44, label = "bias", size = 6) +
    annotate("text", x = 1.075, y = 1.56, label = "higher bias", size = 6) +
    annotate("text", x = 1.075, y = 1.44, label = "with additive", size = 6) +
    geom_segment(aes(x = .9, y = 1.275, xend = 1.4, yend = 1.275),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    theme(aspect.ratio=1, legend.title=element_text(size=8.5)) 
  
  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = sigma_ratio, fill = var_diff_e0)) +
    geom_tile() +
    ylab("Ratio of Environmental Standard Errors") +
    xlab("Difference in Environment Specific Effects (au)") +
    ggtitle("Variance") +
    scale_fill_gradient2(limits = c(-.75, 1.125),
                         breaks = seq(-.75, 1.125, .375),
                         labels = seq(-.75, 1.125, .375)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Difference") +
    annotate("text", x = .75, y = 2.06, label = "lower variance", size = 6) +
    annotate("text", x = .75, y = 1.94, label = "with GxE", size = 6) +
    annotate("text", x = .75, y = .98,
             label = "higher variance", size = 6) +
    annotate("text", x = .75, y = .86,
             label = "with GxE", size = 6) +
    geom_segment(aes(x = .17, y = 1.1, xend = .17, yend = .8),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    geom_segment(aes(x = .17, y = 1.75, xend = .17, yend = 2.1),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    theme(aspect.ratio=1, legend.title=element_text(size=8.5)) 
    
  
  
  return(
    list(
      mse_plot = g1,
      bias_plot = g2,
      variance_plot = g3
    )
  )
  
}

# # first time simulation
# plot1 <- generate_mse_tradeoff_single_site_se_ratio_plot(write_rds = T)
# 
# # first time simulation
plot2 <- generate_mse_tradeoff_single_site_plot(write_rds = T, sigma_ratio = 1)

library(ggpubr)
library(grid)

fig <- ggarrange(
  plot2$bias_plot + rremove("ylab") + rremove("xlab"),
  plot2$variance_plot + rremove("ylab") + rremove("xlab"),
  plot2$mse_plot + rremove("ylab") + rremove("xlab"),
  nrow = 1,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B", "C"),
  vjust = .1
)

#png(file = "~/Documents/paper_gwas_amplification/images/fig1.png", width = 4000, height = 1350, res = 300)
annotate_figure(fig,
                 left = textGrob("Context Specific Estimation Noise (a.u.)", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 14)),
                 bottom = textGrob("Difference in Context Specific Effects (a.u.)", gp = gpar(cex = 1, fontsize = 14), hjust = .575),
                 top = textGrob("MSE Difference in the Case of Context Specific Noise", gp = gpar(cex = 1, fontsize = 15), hjust = .625, vjust = 2)
                )
dev.off()

# fig2 <- ggarrange(
#   plot1$bias_plot + rremove("ylab") + rremove("xlab"),
#   plot1$variance_plot + rremove("ylab") + rremove("xlab"),
#   plot1$mse_plot + rremove("ylab") + rremove("xlab"),
#   nrow = 1,
#   common.legend = TRUE,
#   legend = "right",
#   labels = c("A", "B", "C"),
#   widths = c(1, 1, 1)
# )
# 
# png(file = "~/Documents/paper_gwas_amplification/images/fig2.png", width = 3750, height = 1250, res = 300)
# annotate_figure(fig2,
#                 left = textGrob("Estimation Noise Ratio", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 15)),
#                 bottom = textGrob("Difference in Context Specific Effects (a.u.)", gp = gpar(cex = 1, fontsize = 15), hjust = .55))
# dev.off()

mse_multi_ratio_plot <- generate_multi_ratio_plot()

library(ggpubr)

fig <- ggarrange(
  mse_multi_ratio_plot$mse_plot_rtwo_thirds + rremove("ylab") + rremove("xlab"),
  NULL,
  mse_multi_ratio_plot$mse_plot_r1 + rremove("ylab") + rremove("xlab"),
  NULL,
  mse_multi_ratio_plot$mse_plot_rone_five + rremove("ylab") + rremove("xlab"),
  nrow = 1,
  common.legend = TRUE,
  labels = c("A", "", "B", "", "C"),
  widths = c(1, 0, 1, 0, 1),
  legend = "right",
  vjust = .1
)

#png(file = "~/Documents/paper_gwas_amplification/images/fig2.png", width = 4000, height = 1250, res = 300)
annotate_figure(fig,
                left = textGrob("Estimation Noise in Context A (a.u.)", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 14)),
                bottom = textGrob("Difference in Context Specific Effects (a.u.)", gp = gpar(cex = 1, fontsize = 14), hjust = .625),
                top = textGrob("MSE Difference in the Case of Context Specific Noise", gp = gpar(cex = 1, fontsize = 15), hjust = .625, vjust = 2)
              )
#dev.off()
