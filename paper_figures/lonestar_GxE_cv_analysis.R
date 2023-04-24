'%>%' <- magrittr::'%>%'
# First, load in the data
# I will need to change this for lonestar
#data_dir <- "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results"
data_dir <- "/corral-repl/utexas/Harpak-Lab-GWAS/gxe_decision_rule/final_files/"

create_cv_df <- function(data_folder, trait, n_folds = 5) {
  
  cv_df <- data.frame()
  
  for (fold in 1:n_folds) {
    
    # Training Data
    for (sex in c("male", "female", "all")) {
      
      fold_df <- readr::read_tsv(
        file = glue::glue("{data_folder}/{trait}/{trait}.{sex}_train.1.merged.txt"),
        col_select = c("ID", "BETA", "SE")
      ) %>%
        setNames(c("id", glue::glue("b_train_{sex}_{fold}"), glue::glue("se_train_{sex}_{fold}"))) 
      
      if (nrow(cv_df) == 0) {
        
        cv_df <- fold_df
        
      } else {
        
        cv_df <- cv_df %>% dplyr::inner_join(fold_df, by = "id")
        
      }
      
    }
    
    # Test Data
    for (sex in c("male", "female")) {
      
      fold_df <- readr::read_tsv(
        file = glue::glue("{data_folder}/{trait}/{trait}.{sex}_test.1.merged.txt"),
        col_select = c("ID", "BETA", "SE")
      ) %>%
        setNames(c("id", glue::glue("b_test_{sex}_{fold}"), glue::glue("se_test_{sex}_{fold}"))) 
      
      cv_df <- cv_df %>% dplyr::inner_join(fold_df, by = "id")
      
    }
    
  }
  
  return(cv_df)
  
}

cv_df <- create_cv_df(
  data_dir, "diastolicBP_auto", 1
)

cv_df <- cv_df %>% dplyr::distinct(id, .keep_all = TRUE)

# I think here I should perhaps threshold for true difference in sex effects
cv_df <- cv_df %>%
  dplyr::mutate(
    beta_diff = b_train_male_1 - b_train_female_1,
    se_beta_diff = sqrt(se_train_male_1 ^ 2 + se_train_female_1 ^ 2)
  )

# cv_df <- cv_df %>%
#   dplyr::mutate(
#     bias_corrected_beta_diff = pmax(0, abs(beta_diff) - sqrt(2 / pi) * se_beta_diff)
#   )

cv_df <- cv_df %>%
  tidyr::drop_na()

ash_train_df <- cv_df %>% dplyr::sample_n(250000)

ash_model <- ashr::ash(
  betahat = ash_train_df$beta_diff,
  sebetahat = ash_train_df$se_beta_diff
)

ash_regressed_model_all_snps <- ashr::ash(
  betahat = cv_df$beta_diff,
  sebetahat = cv_df$se_beta_diff,
  fixg = TRUE,
  g = ashr::get_fitted_g(ash_model)
)

cv_df <- cv_df %>%
  dplyr::mutate(
    expected_add_male_bias = -.5 * ashr::get_pm(ash_regressed_model_all_snps),
    expected_add_female_bias = .5 * ashr::get_pm(ash_regressed_model_all_snps),
    expected_add_var = se_train_all_1 ^ 2,
    expected_GxE_male_mse = se_train_male_1 ^ 2,
    expected_GxE_female_mse = se_train_female_1 ^ 2
  )

cv_df <- cv_df %>%
  dplyr::mutate(
    expected_add_male_mse = expected_add_male_bias ^ 2 + expected_add_var,
    expected_add_female_mse = expected_add_female_bias ^ 2 + expected_add_var
  )

cv_df <- cv_df %>%
  dplyr::mutate(
    expected_male_add_mse_advantage = - (expected_add_male_mse - expected_GxE_male_mse),
    expected_female_add_mse_advantage = - (expected_add_female_mse - expected_GxE_female_mse)
  )

cv_df <- cv_df %>%
  dplyr::mutate(
    actual_add_male_mse = (b_train_all_1 - b_test_male_1) ^ 2,
    actual_add_female_mse = (b_train_all_1 - b_test_female_1) ^ 2,
    actual_GxE_male_mse = (b_train_male_1 - b_test_male_1) ^ 2,
    actual_GxE_female_mse = (b_train_female_1 - b_test_female_1) ^ 2
  )

cv_df <- cv_df %>%
  dplyr::mutate(
    actual_male_add_mse_advantage = - (actual_add_male_mse - actual_GxE_male_mse),
    actual_female_add_mse_advantage = - (actual_add_female_mse - actual_GxE_female_mse)
  )

readr::write_csv(cv_df, "cv_diastolicBP_fold1_ash_full.csv")

# samp_df <- cv_df %>%
#   dplyr::sample_n(100000)

# at this point I have to divide up the data into discrete points myself
# I think that this should be done locally
