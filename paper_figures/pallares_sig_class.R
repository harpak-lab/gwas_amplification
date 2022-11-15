library(ggplot2)
'%>%' <- magrittr::'%>%'

summary_table <- read.delim('data/SummaryTable_allsites_12Nov20.txt')

# replace 0 p-values with small numbers
summary_table <- summary_table %>%
  dplyr::select(c(site, pval_CTRL, pval_HS, coef_CTRL, coef_HS, sig_cat)) %>%
  dplyr::mutate(
    pval_CTRL = pmax(1e-12, pval_CTRL),
    pval_HS = pmax(1e-12, pval_HS)
  )

# construct std error estimates from coefficients and p-values
summary_table <- summary_table %>%
  dplyr::mutate(
    std_error_ctrl = abs(coef_CTRL) / qnorm((2 - pval_CTRL) / 2),
    std_error_hs = abs(coef_HS) / qnorm((2 - pval_HS) / 2)
  )

summary_table_signif <- summary_table %>%
  dplyr::filter(sig_cat != "NS")

summary_table_signif <- summary_table_signif %>%
  dplyr::mutate(sig_diff = abs(coef_CTRL - coef_HS))

plot_line <- gamp::get_sigma_ratio_line(1.2)

# Rename here
plot_colors <- c("green", "orange", "blue")
names(plot_colors) <- c(
  "CTRL", "HS", "shared"
)

ggplot(data = summary_table_signif) +
  geom_point(size = .5, aes(x = sig_diff, y = std_error_hs, color = sig_cat), alpha = .4) +
  geom_segment(
    aes(x = 0.075, y = .075 * plot_line$slope, xend = .325, yend = .325 * plot_line$slope),
    linetype = "dashed"
  ) +
  ylab("High Sugar Effect Standard Error") +
  xlab("Estimated Difference in High Sugar and Control Effects") +
  ggtitle("MSE Decision Boundary for Fly Data") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color="Pallares Classification") +
  scale_color_manual(values = plot_colors)

