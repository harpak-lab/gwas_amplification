# Here, I'd like to examine Carrie's cross validation results
# I wan't to see what I can find from them / think about how cross-validation might work
'%>%' <- magrittr::'%>%'

create_cv_df <- function(data_folder, trait, n_folds = 5) {
  
  cv_df <- data.frame()
  
  for (fold in 1:n_folds) {
    
    for (sex in c("male", "female", "both_sex")) {
      
      if (fold == 1 && sex == "male") {
        
        fold_df <- readr::read_tsv(
          file = glue::glue("{data_folder}/{trait}/{trait}_PGS{fold}/{sex}_train.{trait}.glm.linear"),
          col_select = c("ID", "BETA", "SE"),
          n_max = 2000000
        ) %>%
          setNames(c("id", glue::glue("b_{sex}_{fold}"), glue::glue("se_{sex}_{fold}"))) 
        
      } else {
        
        fold_df <- readr::read_tsv(
          file = glue::glue("{data_folder}/{trait}/{trait}_PGS{fold}/{sex}_train.{trait}.glm.linear"),
          col_select = c("BETA", "SE"),
          n_max = 2000000
        ) %>%
          setNames(c(glue::glue("b_{sex}_{fold}"), glue::glue("se_{sex}_{fold}"))) 
        
      }
      

      
      if (nrow(cv_df) == 0) {
        
        cv_df <- fold_df
        
      } else {
        
        cv_df <- cbind(cv_df, fold_df)
        
      }

    }
    
  }
  
  return(cv_df)
  
}

get_sum_stats_loo <- function(cv_df, left_out_fold, n_folds = 5) {
  
  loo_df <- data.frame(id = cv_df$id)
  n_folds_used <- n_folds - 1
  
  for (sex in c("male", "female", "both_sex")) {
    
    sex_b_fold_names <- stringr::str_c(glue::glue("b_{sex}_"), seq(1, 5)[-left_out_fold])
    loo_b <- cv_df %>%
      dplyr::select(sex_b_fold_names) %>%
      as.matrix() %>%
      rowMeans()
    
    sex_se_fold_names <- stringr::str_c(glue::glue("se_{sex}_"), seq(1, 5)[-left_out_fold])
    loo_se_mat <- cv_df %>%
      dplyr::select(sex_se_fold_names) %>%
      as.matrix()
    
    loo_se <- sqrt((loo_se_mat ^ 2) %*% rep((1 / n_folds_used) ^ 2, n_folds_used))
    
    loo_df[[glue::glue("b_{sex}")]] <- loo_b
    loo_df[[glue::glue("se_{sex}")]] <- loo_se
    
  }
  
  return(loo_df)
  
}

get_dist_from_line <- function(slope, intercept, x, y) {
  
  line_y <- intercept + slope * x
  dist <- y - line_y
  return(dist)
  
}

cv_df <- create_cv_df(data_folder = "~/Downloads", trait = "height")

# Now, I can iterate through each fold, evaluate the decision rule, and then record stats
n_folds <- 5
se_ratio <- mean(cv_df$se_male_1 / cv_df$se_female_1)
dec_rule <- gamp:::get_sigma_ratio_line(se_ratio)
# we can focus on predicting in females

errors_df <- data.frame()

for (i in 1:n_folds) {
  
  print(glue::glue("Fold {i}"))
  
  # get sumstats from other folds
  sum_stat_loo <- get_sum_stats_loo(cv_df, left_out_fold = i)
  t_stats <- abs((sum_stat_loo$b_male - sum_stat_loo$b_female) / sqrt(
    sum_stat_loo$se_male ^ 2 + sum_stat_loo$se_female ^ 2)
  )
  
  t_stats_p <- dnorm(t_stats)
  fx_diff <- dplyr::if_else(
    t_stats_p < .05, abs(sum_stat_loo$b_male - sum_stat_loo$b_female), 0
  )
  
  # an interesting filter we could perhaps add here is a t-test for difference in means
  # I'm not sure if this would help
  distances_t <- get_dist_from_line(
    slope = dec_rule$slope, 
    intercept = dec_rule$intercept,
    x = fx_diff, 
    y = sum_stat_loo$se_female
  )
  
  distances <- get_dist_from_line(
    slope = dec_rule$slope, 
    intercept = dec_rule$intercept,
    x = abs(sum_stat_loo$b_male - sum_stat_loo$b_female), 
    y = sum_stat_loo$se_female
  )
  
  # Now, get error for predicting females from additive and from female specific
  additive_error <- abs(cv_df[[glue::glue("b_female_{i}")]] - sum_stat_loo$b_both_sex)
  GxE_error <- abs(cv_df[[glue::glue("b_female_{i}")]] - sum_stat_loo$b_female)
  
  fold_error_df <- data.frame(
    distance = distances,
    distance_t = distances_t,
    additive_error = additive_error,
    GxE_error = GxE_error
  )
  
  errors_df <- rbind(errors_df, fold_error_df)
  
}

center_errors_df <- errors_df %>%
  dplyr::filter(
    distance > quantile(distance, .01) & distance < quantile(distance, .99)
  )

center_errors_df <- center_errors_df %>%
  dplyr::mutate(abs_error_diff = GxE_error - additive_error)

center_errors_df <- center_errors_df %>%
  dplyr::mutate(decision_rule = dplyr::if_else(distance >= 0, "additive", "GxE"))

summ_center_errors_df <- center_errors_df %>%
  dplyr::group_by(decision_rule) %>%
  dplyr::summarize(
    total = dplyr::n(),
    additive_mse = round(mean(additive_error ^ 2), 5),
    GxE_mse = round(mean(GxE_error ^ 2), 5)
  )


