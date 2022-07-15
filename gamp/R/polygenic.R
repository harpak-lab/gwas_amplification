#' Get Observation Noise Based on Heritability
#'
#' @param maf_vec vector containing the minor allele frequency of each snp
#' @param fx_vec vector containing true effect size in each snp
#' @param h2 number indiciating desired heritability
#'
#' @return environmental variance
#'
get_theo_noise_from_h2 <- function(
  maf_vec, fx_vec, h2
) {

  var_G <- sum((fx_vec ^ 2) * (2 * (1 - maf_vec) * maf_vec))
  var_E <- (var_G - h2 * var_G) / h2
  return(sqrt(var_E))

}


#' Simulate polygenic data from two environments
#'
#' The model assumes that some outcome y is generated as a linear function
#' of the true effects in each environment
#'
#' @param fx_mat n x 2m matrix of true effects
#' @param n_e0,n_e1 number of individuals in each environment
#' @param sigma_e0,sigma_e1 observation noise in each environment
#'
#' @return dataframe of simulated genotypes and effects
#'
simulate_polygenic_data_2env <- function(
  fx_mat,
  n_e0,
  n_e1,
  e0_h2,
  e1_h2,
  maf_simulator
) {

  n_snps <- nrow(fx_mat)
  sim_df <- data.frame(e = c(rep("e0", n_e0), rep("e1", n_e1)))

  # simulate minor allele frequencies all at once for efficiency
  maf_vec <- do.call(maf_simulator, list(n_snps))
  y_e0 <- numeric(n_e0)
  y_e1 <- numeric(n_e1)

  for (i in 1:n_snps) {

    genotypes <- rbinom(n = n_e0 + n_e1, size = 2, prob = maf_vec[i])
    y_e0 <- y_e0 + genotypes[1:n_e0] * fx_mat[i, 1]
    y_e1 <- y_e1 + genotypes[(n_e0 + 1):(n_e0 + n_e1)] * fx_mat[i, 2]
    sim_df[[glue::glue('snp_{i}_g')]] <- genotypes

  }

  sigma_e0 <- get_theo_noise_from_h2(
    maf_vec, fx_mat[, 1], e0_h2
  )
  sigma_e1 <- get_theo_noise_from_h2(
    maf_vec, fx_mat[, 2], e1_h2
  )

  y_e0 <- y_e0 + rnorm(n = n_e0, sd = sigma_e0)
  y_e1 <- y_e1 + rnorm(n_e1, sd = sigma_e1)

  sim_df$y <- c(y_e0, y_e1)

  return(sim_df)

}

#' Get predictions from an additive model
#'
#' @param train_df dataframe to fit the additive model
#' @param test_df dataframe on which to make polygenic outcome predictions
#'
#' @return dataframe of predictions and effect size estimates
#'
get_test_preds_additive <- function(
  train_df,
  test_df
) {

  # 1 column for y and 1 column for e
  n_snps <- ncol(train_df) - 2
  n_pred <- nrow(test_df)

  pred_y <- numeric(n_pred)
  est_coefs <- numeric(n_snps)

  for (i in 1:n_snps) {

    snp_name <- glue::glue("snp_{i}_g")
    train_snp_genotypes <- train_df[[snp_name]]
    snp_lm <- lm(train_df$y ~ train_snp_genotypes)
    snp_coef <- coef(summary(snp_lm))['train_snp_genotypes', 'Estimate']

    pred_y <- pred_y + snp_coef * test_df[[snp_name]]
    est_coefs[i] <- snp_coef

  }

  return(list(
    pred = pred_y,
    est = est_coefs
  ))

}

#' Get predictions from a GxE model
#'
#' @param train_df dataframe to fit the additive model
#' @param test_df dataframe on which to make polygenic outcome predictions
#'
#' @return dataframe of predictions and effect size estimates in both environments
#'
get_test_preds_GxE <- function(
  train_df,
  test_df
) {

  # 1 column for y and 1 column for e
  n_snps <- ncol(train_df) - 2

  train_df_e0 <- train_df %>%
    dplyr::filter(e == "e0")

  train_df_e1 <- train_df %>%
    dplyr::filter(e == "e1")

  test_df_e0 <- test_df %>%
    dplyr::filter(e == "e0")

  test_df_e1 <- test_df %>%
    dplyr::filter(e == "e1")

  n_pred_e0 <- nrow(test_df_e0)
  n_pred_e1 <- nrow(test_df) - n_pred_e0

  pred_y_e0 <- numeric(n_pred_e0)
  pred_y_e1 <- numeric(n_pred_e1)

  est_coefs_e0 <- numeric(n_snps)
  est_coefs_e1 <- numeric(n_snps)

  for (i in 1:n_snps) {

    snp_name <- glue::glue("snp_{i}_g")

    train_snp_genotypes_e0 <- train_df_e0[[snp_name]]
    snp_e0_lm <- lm(train_df_e0$y ~ train_snp_genotypes_e0)
    snp_e0_coef <- coef(summary(snp_e0_lm))['train_snp_genotypes_e0', 'Estimate']
    pred_y_e0 <- pred_y_e0 + snp_e0_coef * test_df_e0[[snp_name]]

    train_snp_genotypes_e1 <- train_df_e1[[snp_name]]
    snp_e1_lm <- lm(train_df_e1$y ~ train_snp_genotypes_e1)
    snp_e1_coef <- coef(summary(snp_e1_lm))['train_snp_genotypes_e1', 'Estimate']
    pred_y_e1 <- pred_y_e1 + snp_e1_coef * test_df_e1[[snp_name]]

    est_coefs_e0[i] <- snp_e0_coef
    est_coefs_e1[i] <- snp_e1_coef

  }

  return(list(
    preds = c(pred_y_e0, pred_y_e1),
    est_e0 = est_coefs_e0,
    est_e1 = est_coefs_e1
  ))

}

#' Get predictions from a mash model
#'
#' @param train_df dataframe to fit the additive model
#' @param test_df dataframe on which to make polygenic outcome predictions
#' @param corr_range vector indicating range of correlation for covariance matrices.
#' Must be beteween 0 and 1.
#' @param amp_range vector indicating range of amplification of covariance matrices.
#'
#' @return dataframe of predictions and effect size estimates in both environments
#'
get_test_preds_mash <- function(
  train_df,
  test_df,
  corr_range = seq(from = -1, to = 1, by = (1/3)),
  amp_range = c(1.5, 2)
) {

  # 1 column for y and 1 column for e
  n_snps <- ncol(train_df) - 2

  # mash setup
  Bhat <- matrix(nrow = n_snps, ncol = 2)
  Shat <- matrix(nrow = n_snps, ncol = 2)

  train_df_e0 <- train_df %>%
    dplyr::filter(e == "e0")

  train_df_e1 <- train_df %>%
    dplyr::filter(e == "e1")

  test_df_e0 <- test_df %>%
    dplyr::filter(e == "e0")

  test_df_e1 <- test_df %>%
    dplyr::filter(e == "e1")

  for (i in 1:n_snps) {

    snp_name <- glue::glue("snp_{i}_g")

    train_snp_genotypes_e0 <- train_df_e0[[snp_name]]
    snp_e0_lm <- lm(train_df_e0$y ~ train_snp_genotypes_e0)
    snp_e0_coef <- coef(summary(snp_e0_lm))['train_snp_genotypes_e0', 'Estimate']
    snp_e0_se <- coef(summary(snp_e0_lm))['train_snp_genotypes_e0', 'Std. Error']

    train_snp_genotypes_e1 <- train_df_e1[[snp_name]]
    snp_e1_lm <- lm(train_df_e1$y ~ train_snp_genotypes_e1)
    snp_e1_coef <- coef(summary(snp_e1_lm))['train_snp_genotypes_e1', 'Estimate']
    snp_e1_se <- coef(summary(snp_e1_lm))['train_snp_genotypes_e1', 'Std. Error']

    Bhat[i, 1] <- snp_e0_coef
    Bhat[i, 2] <- snp_e1_coef

    Shat[i, 1] <- snp_e0_se
    Shat[i, 2] <- snp_e1_se

  }

  mash_fit <- fit_mash_gwas_2env(
    Bhat = Bhat,
    Shat = Shat,
    corr_range = corr_range,
    amp_range = amp_range
  )

  mash_est_fx <- ashr::get_pm(mash_fit)

  n_pred_e0 <- nrow(test_df_e0)
  n_pred_e1 <- nrow(test_df) - n_pred_e0

  pred_y_e0 <- numeric(n_pred_e0)
  pred_y_e1 <- numeric(n_pred_e1)

  for (i in 1:n_snps) {

    snp_name <- glue::glue("snp_{i}_g")
    pred_y_e0 <- pred_y_e0 + mash_est_fx[i, 1] * test_df_e0[[snp_name]]
    pred_y_e1 <- pred_y_e1 + mash_est_fx[i, 2] * test_df_e1[[snp_name]]

  }

  return(list(
    preds = c(pred_y_e0, pred_y_e1),
    est_e0 = mash_est_fx[,1],
    est_e1 = mash_est_fx[,2]
  ))

}

#' Simulate correlation between prediction of an outcome and polygenic
#' scores for various models
#'
#' @param n_e0,n_e1 number of individuals in each environment
#' @param n_snps number of SNPs
#' @param maf_simulator function that takes number of snps as an argument and
#' returns minor allele frequencies for each snp
#' @param sigma_e0,sigma_e1 additive noise standard deviation in each environment
#' @param Sigma list of covariance matrices from which true effects are drawn
#' @param pi vector of probabilities corresponding to Sigma
#' @param num_sims number of simulations to run
#' @param s (optional) level of sparsity of true effects. Defaults to 0
#'
#' @return
#' @export
#'
#' @examples
simulate_pgs_corr <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0
) {

  mash_corr_vec <- numeric(num_sims)
  additive_corr_vec <- numeric(num_sims)
  GxE_corr_vec <- numeric(num_sims)

  mash_mse_vec <- numeric(num_sims)
  additive_mse_vec <- numeric(num_sims)
  GxE_mse_vec <- numeric(num_sims)

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    # obtain polygenic data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0,
      n_e1,
      e0_h2,
      e1_h2,
      maf_simulator
    )

    ## 50% of the sample size
    smp_size <- floor(0.5 * (n_e0 + n_e1))

    train_ind <- sample(seq_len(n_e0 + n_e1), size = smp_size)

    train_df <- polygenic_df[train_ind, ]
    test_df <- polygenic_df[-train_ind, ]

    additive_preds <- get_test_preds_additive(train_df, test_df)
    GxE_preds <- get_test_preds_GxE(train_df, test_df)
    mash_preds <- get_test_preds_mash(train_df, test_df)

    # here, I can change the error metric to MSE
    additive_corr_vec[i] <- cor(additive_preds$pred, test_df$y)
    GxE_corr_vec[i] <- cor(GxE_preds$pred, test_df$y)
    mash_corr_vec[i] <- cor(mash_preds$pred, test_df$y)

    additive_mse_vec[i] <- mean(
      (fx_mat - matrix(data = c(additive_preds$est, additive_preds$est), ncol = 2)) ^ 2
    )
    GxE_mse_vec[i] <- mean(
      (fx_mat - matrix(data = c(GxE_preds$est_e0, GxE_preds$est_e1), ncol = 2)) ^ 2
    )
    mash_mse_vec[i] <- mean(
      (fx_mat - matrix(data = c(mash_preds$est_e0, mash_preds$est_e1), ncol = 2)) ^ 2
    )

  }

  results_list <- list(
    additive_cor = mean(additive_corr_vec),
    GxE_cor = mean(GxE_corr_vec),
    mash_cor = mean(mash_corr_vec),
    additive_mse = mean(additive_mse_vec),
    GxE_mse = mean(GxE_mse_vec),
    mash_mse = mean(mash_mse_vec)
  )

  return(results_list)

}
