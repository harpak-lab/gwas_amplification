tictoc::tic()
sim_beta_diff_list <- gamp:::sim_beta_hat_diff(
  diff_range = seq(from = -1, to = 1, length.out = 40), 
  sd_range = seq(from = .01, to = 1, length.out = 40), 
  pts_per_sim = 1250
)
tictoc::toc()

ggplot(data = sim_beta_diff_list$beta_hat_diff_df, aes(x = mean_abs_diff, y = sd, fill = e1_add_mse - e1_GxE_mse)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  ylab("Environment Specific Standard Error") +
  xlab("Estimated Difference in Environment Specific Effects") +
  ggtitle("MSE") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_text(size=7))

ggplot(data = sim_beta_diff_list$beta_hat_bc_sum_sim_df, aes(x = mean_abs_bc_diff, y = sd, fill = e1_add_mse - e1_GxE_mse)) +
  geom_tile() +
  scale_fill_gradient2() +
  coord_fixed() +
  xlab("Bias Corrected Difference in Effect Size") +
  ylab("Environment Specific Standard Error") +
  ggtitle("MSE") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_text(size=7))

