'%>%' <- magrittr::'%>%'

df_list <- list()

num_dfs <- 1

for (i in seq(0, 11022, 1002)) {
  
  df_list[[num_dfs]] <- readr::read_csv(
    glue::glue("~/Documents/gwas_amplification/paper_figures/pallares_cv_results/pallares_cv_start_{i+1}_stop_{i+1002}.csv")
  )
  num_dfs <- num_dfs + 1
  
}

pallares_cv_df <- dplyr::bind_rows(df_list, .id = "column_label") %>%
  dplyr::select(-column_label)

# Now, I have to come up with a method for the decision rule
# I think that I should just compare the bias and the variance between the two methods
# and this should work

library(dplyr)

# want to come back here after I try this for the GxSex data
# I'm curious what would happen if I threshold for any difference in envs

pallares_cv_df <- pallares_cv_df %>%
  dplyr::mutate(
    env_beta_diff = b_c_train - b_hs_train,
    env_se_diff = sqrt(se_c_train ^ 2 + se_hs_train ^ 2)
    #env_diff_z_stat = (b_c_train - b_hs_train) / (se_c_train ^ 2 + se_hs_train ^ 2)
  )

pallares_cv_df <- pallares_cv_df %>% tidyr::drop_na()

# now, want to train an ash model
pallares_ash_model <- ashr::ash(
  betahat = pallares_cv_df$env_beta_diff, 
  sebetahat = pallares_cv_df$env_se_diff,
  mixcompdist = "normal",
  method = "shrink"
)

pallares_cv_df <- pallares_cv_df %>%
  dplyr::mutate(
    env_diff_pm = ashr::get_pm(pallares_ash_model)
  )

#pval_thresh <- .05

pallares_cv_df <- pallares_cv_df %>%
  mutate(
    expected_hs_mse = se_hs_train ^ 2,
    expected_c_mse = se_c_train ^ 2,
    expected_t_hs_mse = se_t_train ^ 2 + (env_diff_pm / 2) ^ 2,
    expected_t_c_mse = se_t_train ^ 2 + (env_diff_pm / 2) ^ 2
  )

pallares_cv_df <- pallares_cv_df %>%
  mutate(
    actual_hs_mse = (b_hs_train - b_hs_test) ^ 2,
    actual_c_mse = (b_c_train - b_c_test) ^ 2,
    actual_t_hs_mse = (b_t_train - b_hs_test) ^ 2,
    actual_t_c_mse = (b_t_train - b_c_test) ^ 2
  )

pallares_cv_df <- pallares_cv_df %>%
  mutate(
    expected_c_additive_benefit = -(expected_t_c_mse - expected_c_mse),
    expected_hs_additive_benefit = -(expected_t_hs_mse - expected_hs_mse),
    actual_c_additive_benefit = -(actual_t_c_mse - actual_c_mse),
    actual_hs_additive_benefit = -(actual_t_hs_mse - actual_hs_mse)
  )

# Exclusion of outliers
pallares_cv_df <- pallares_cv_df %>%
  dplyr::filter(abs(se_c_test) < 2 & abs(se_hs_test) < 2)

pallares_cv_df$decile <- cut(
  pallares_cv_df$expected_c_additive_benefit,
  quantile(pallares_cv_df$expected_c_additive_benefit, probs = seq(0, 1, .1)),
  include.lowest=TRUE,
  labels=FALSE
)

c_summ_df <- pallares_cv_df %>%
  dplyr::group_by(decile) %>%
  dplyr::summarize(
    weighted_dec_advantage_add_c = weighted.mean(actual_c_additive_benefit, 1 / (se_c_test ^ 2)),
    unweighted_dec_advantage_add_c = mean(actual_c_additive_benefit)
  )
