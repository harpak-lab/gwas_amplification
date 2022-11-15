# want to make a plot with the Pallares data
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

summary_table <- summary_table %>%
  dplyr::filter(sig_cat != "NS")

summary_table <- summary_table %>%
  dplyr::mutate(sig_diff = abs(coef_CTRL - coef_HS))

ggplot(data = summary_table, aes(x = abs(coef_HS), y = std_error_hs)) +
  geom_point(alpha = .1, color = "blue") +
  geom_segment(aes(x = .1, y = .0678, xend = .6, yend = 0.4068), linetype = "dashed")

