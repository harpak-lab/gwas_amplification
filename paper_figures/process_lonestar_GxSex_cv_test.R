# first, I want to read in the data from lonestar
# then I can make a grid
'%>%' <- magrittr::'%>%'
cv_df <- readr::read_csv(
  file = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results/cv_testosterone_fold1_ash_full.csv"
)

dec_rule <- gamp:::get_sigma_ratio_slope(
  r = median(cv_df$se_train_male_1 / cv_df$se_train_female_1), n = 147477, m = 143937
)


male_cv_df <- cv_df %>%
  dplyr::select(se_train_male_1, expected_add_male_bias, actual_male_add_mse_advantage) %>%
  dplyr::mutate(fx_diff = 2 * abs(expected_add_male_bias)) %>%
  dplyr::rename(se = se_train_male_1, mse_diff = actual_male_add_mse_advantage) %>%
  dplyr::select(-expected_add_male_bias)

rm("cv_df")
gc()

male_cv_df_samp <- male_cv_df %>%
  dplyr::sample_n(100000)

library(ggplot2)

ggplot(data = male_cv_df_samp, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .025, size = 1, color = "orange") +
  ylab("Male Standard Error (ng/dl)") +
  xlab("Difference in Environment Specific Effects (ng/dl)") +
  ggtitle("Testosterone MSE") +
  geom_segment(
    aes(x = 0, y = 0, xend = .1, yend = .1 * dec_rule), linetype = "dashed"
  ) +
  labs(fill=latex2exp::TeX("Difference $(ng/dl)^{2}$")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(0, .1) +
  ylim(0, .5)


# the code below is for making the correlation plot
male_cv_df <- cv_df %>%
  dplyr::filter(
    expected_male_add_mse_advantage > quantile(expected_male_add_mse_advantage, probs = .025) &
      expected_male_add_mse_advantage < quantile(expected_male_add_mse_advantage, probs = .975)
  )

# I think it would be useful to perform a boxing strategy here
# and see how the results look
# this will help me understand where the real differences lie

# male_cv_df$decile <- cut(
#   male_cv_df$expected_male_add_mse_advantage,
#   quantile(male_cv_df$expected_male_add_mse_advantage, probs = seq(0, 1, .1)),
#   include.lowest=TRUE,
#   labels=FALSE
# )

male_cv_df <- male_cv_df %>%
  dplyr::mutate(
    percentile = cut(
      expected_male_add_mse_advantage,
      seq(min(expected_male_add_mse_advantage), max(expected_male_add_mse_advantage), length.out = 10),
      include.lowest = TRUE,
      labels = FALSE
    )
  )


male_summ_df <- male_cv_df %>%
  dplyr::group_by(percentile) %>%
  dplyr::summarize(
    weighted_dec_advantage_add_male = weighted.mean(actual_male_add_mse_advantage, 1 / (se_test_male_1 ^ 2)),
    unweighted_dec_advantage_add_male = mean(actual_male_add_mse_advantage),
    weighted_sd_advantage_add_male = sqrt(Hmisc::wtd.var(actual_male_add_mse_advantage, 1 / (se_test_male_1 ^ 2))),
    expected_advantage = mean(expected_male_add_mse_advantage),
    samp_size = dplyr::n(),
    sd = var(actual_male_add_mse_advantage)
  )

library(ggplot2)

ggplot(data = male_summ_df, 
       aes(x = expected_advantage, y = unweighted_dec_advantage_add_male,
           ymin = unweighted_dec_advantage_add_male - sqrt(sd / samp_size),
           ymax = unweighted_dec_advantage_add_male + sqrt(sd / samp_size)
       )) +
  geom_point() +
  geom_errorbar() +
  xlab("Expected Male MSE Advantage") +
  ylab("Actual Male MSE Advantage") +
  ggtitle("Testosterone")

plot(
  male_summ_df$expected_advantage, 
  male_summ_df$weighted_dec_advantage_add_male,
  xlab = "Expected Male MSE Advantage",
  ylab = "Actuale Male MSE Advantage"
)




