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
  maf_vec
) {

  total_n <- n_e0 + n_e1
  n_snps <- nrow(fx_mat)
  sim_df <- data.frame(e = c(rep("e0", n_e0), rep("e1", n_e1)))

  genotypes_mat <- matrix(
    data = rbinom(n = n_snps * total_n, size = 2, prob = maf_vec),
    nrow = total_n,
    ncol = n_snps,
    byrow = TRUE
  )

  sigma_e0 <- get_theo_noise_from_h2(
    maf_vec, fx_mat[, 1], e0_h2
  )
  sigma_e1 <- get_theo_noise_from_h2(
    maf_vec, fx_mat[, 2], e1_h2
  )

  y_e0 <- genotypes_mat[1:n_e0, ] %*% fx_mat[, 1] + rnorm(n_e0, sd = sigma_e0)
  y_e1 <- genotypes_mat[(n_e0 + 1):total_n, ] %*% fx_mat[, 2] + rnorm(n_e1, sd = sigma_e1)

  colnames(genotypes_mat) <- paste0("g", 1:n_snps)
  sim_df <- cbind(sim_df, genotypes_mat)

  sim_df$y <- c(y_e0, y_e1)

  return(sim_df)

}

#' Get predictions from an additive model
#'
#' @param train_df dataframe to fit the additive model
#'
#' @return dataframe of predictions and effect size estimates
#'
estimate_additive_models <- function(
  train_df
) {

  print("Estimating Additive Model")

  # 1 column for y and 1 column for e
  n_snps <- ncol(train_df) - 2

  est_coefs <- numeric(n_snps)
  pvals <- numeric(n_snps)

  design_mat <- matrix(data = 1, nrow = nrow(train_df), ncol = 2)

  snp_names <- paste0("g", 1:n_snps)

  for (i in 1:n_snps) {

    design_mat[, 2] <- train_df[[snp_names[i]]]
    # fast linear regression
    snp_lm <- RcppEigen::fastLm(design_mat, train_df$y)
    pval <- 2 * (1 - pt(
      q = abs(snp_lm$coefficients[2]) / snp_lm$se[2],
      df = snp_lm$df.residual)
    )
    snp_coef <- snp_lm$coefficients[2]
    est_coefs[i] <- snp_coef
    pvals[i] = pval

  }

  return(list(
    pval = pvals,
    est = est_coefs
  ))

}

#' Get Additive Effect Size Estimates
#'
#' @param additive_models list of additive models
#' @param pval_thresh pval threshold
#'
#' @return matrix of estimates
#'
get_ests_additive <- function(additive_models, pval_thresh) {

  estimates <- additive_models$est * (additive_models$pval < pval_thresh)
  est_mat <- matrix(data = c(estimates, estimates), ncol = 2)
  return(est_mat)

}

#' Get predictions from a GxE model
#'
#' @param train_df dataframe to fit the additive model
#'
#' @return dataframe of predictions and effect size estimates in both environments
#'
estimate_GxE_models <- function(
  train_df
) {

  print("Estimating GxE Model")
  # 1 column for y and 1 column for e
  n_snps <- ncol(train_df) - 2

  train_df_e0 <- train_df %>%
    dplyr::filter(e == "e0")

  train_df_e1 <- train_df %>%
    dplyr::filter(e == "e1")

  est_coefs_e0 <- numeric(n_snps)
  est_coefs_e1 <- numeric(n_snps)

  est_sds_e0 <- numeric(n_snps)
  est_sds_e1 <- numeric(n_snps)

  pvals_e0 <- numeric(n_snps)
  pvals_e1 <- numeric(n_snps)

  snp_names <- paste0("g", 1:n_snps)

  design_mat_e0 <- matrix(data = 1, nrow = nrow(train_df_e0), ncol = 2)
  design_mat_e1 <- matrix(data = 1, nrow = nrow(train_df_e1), ncol = 2)

  for (i in 1:n_snps) {

    design_mat_e0[, 2] <- train_df_e0[[snp_names[i]]]
    snp_e0_lm <- RcppEigen::fastLm(design_mat_e0, train_df_e0$y)
    pval_e0 <- 2 * (1 - pt(
      q = abs(snp_e0_lm$coefficients[2]) / snp_e0_lm$se[2],
      df = snp_e0_lm$df.residual)
    )
    snp_e0_coef <- snp_e0_lm$coefficients[2]
    snp_e0_se <- snp_e0_lm$se[2]

    design_mat_e1[, 2] <- train_df_e1[[snp_names[i]]]
    snp_e1_lm <- RcppEigen::fastLm(design_mat_e1, train_df_e1$y)
    pval_e1 <- 2 * (1 - pt(
      q = abs(snp_e1_lm$coefficients[2]) / snp_e1_lm$se[2],
      df = snp_e1_lm$df.residual)
    )
    snp_e1_coef <- snp_e1_lm$coefficients[2]
    snp_e1_se <- snp_e1_lm$se[2]

    est_coefs_e0[i] <- snp_e0_coef
    est_coefs_e1[i] <- snp_e1_coef

    est_sds_e0[i] <- snp_e0_se
    est_sds_e1[i] <- snp_e1_se

    pvals_e0[i] <- pval_e0
    pvals_e1[i] <- pval_e1

  }

  return(list(
    est_e0 = est_coefs_e0,
    sd_e0 = est_sds_e0,
    pval_e0 = pvals_e0,
    est_e1 = est_coefs_e1,
    sd_e1 = est_sds_e1,
    pval_e1 = pvals_e1
  ))

}

#' Get GxE Effect Size Estimates
#'
#' @param gxe_models list of GxE models
#' @param pval_thresh pval threshold
#'
#' @return matrix of estimates
#'
get_ests_GxE <- function(gxe_models, pval_thresh) {

  estimates_e0 <- gxe_models$est_e0 * (gxe_models$pval_e0 < pval_thresh)
  estimates_e1 <- gxe_models$est_e1 * (gxe_models$pval_e1 < pval_thresh)
  est_mat <- matrix(data = c(estimates_e0, estimates_e1), ncol = 2)
  return(est_mat)

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
  Bhat,
  Shat,
  test_df,
  corr_range = seq(from = -1, to = 1, by = (1/3)),
  amp_range = c(1.5, 2),
  fitted_g = NULL
) {

  # 1 column for y and 1 column for e
  n_snps <- nrow(Bhat)
  snp_names <- paste0("g", 1:n_snps)

  test_mat_e0 <- test_df %>%
    dplyr::filter(e == "e0") %>%
    dplyr::select(snp_names) %>%
    as.matrix()

  test_mat_e1 <- test_df %>%
    dplyr::filter(e == "e1") %>%
    dplyr::select(snp_names) %>%
    as.matrix()

  mash_fit <- fit_mash_gwas_2env(
    Bhat = Bhat,
    Shat = Shat,
    corr_range = corr_range,
    amp_range = amp_range,
    fitted_g = fitted_g
  )

  mash_est_fx <- ashr::get_pm(mash_fit)
  mash_lfsr <- ashr::get_lfsr(mash_fit)

  pred_y_e0 <- test_mat_e0 %*% mash_est_fx[, 1]
  pred_y_e1 <- test_mat_e1 %*% mash_est_fx[, 2]

  return(list(
    preds = c(pred_y_e0, pred_y_e1),
    est_e0 = mash_est_fx[,1],
    est_e1 = mash_est_fx[,2],
    lfsr_e0 = mash_lfsr[, 1],
    lfsr_e1 = mash_lfsr[, 2]
  ))

}

mse <- function(Y, X) {

  mean((Y - X) ^ 2)

}

evaluate_selection <- function(fx_mat, est_mat, sel_mat) {


  est_mat[!sel_mat] <- 0

  true_signals <- fx_mat != 0
  est_signals <- est_mat != 0

  power <- sum(true_signals & est_signals) / sum(true_signals)
  type_1_error <- sum(!true_signals & est_signals) / sum(!true_signals)
  accuracy <- sum(est_mat == true_signals) / (2 * nrow(true_signals))

  if (sum(sel_mat) > 0) {

    mse <- mean((est_mat[sel_mat] - fx_mat[sel_mat]) ^ 2)

  } else {
    mse <- 0
  }


  return(
    list(
      power = power,
      type_1_error = type_1_error,
      accuracy = accuracy,
      mse = mse
    )
  )

}

protected_cor <- function(x, y) {

  if (all(x == 0)) {

    return(0)

  } else {

    return(cor(x, y))

  }

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
  mash_mse_vec <- numeric(num_sims)

  additive_corr_vec_pv1 <- numeric(num_sims)
  GxE_corr_vec_pv1 <- numeric(num_sims)

  additive_mse_vec_pv1 <- numeric(num_sims)
  GxE_mse_vec_pv1 <- numeric(num_sims)

  additive_corr_vec_pv2 <- numeric(num_sims)
  GxE_corr_vec_pv2 <- numeric(num_sims)

  additive_mse_vec_pv2 <- numeric(num_sims)
  GxE_mse_vec_pv2 <- numeric(num_sims)

  additive_corr_vec_pv3 <- numeric(num_sims)
  GxE_corr_vec_pv3 <- numeric(num_sims)

  additive_mse_vec_pv3 <- numeric(num_sims)
  GxE_mse_vec_pv3 <- numeric(num_sims)

  additive_corr_vec_pv4 <- numeric(num_sims)
  GxE_corr_vec_pv4 <- numeric(num_sims)

  additive_mse_vec_pv4 <- numeric(num_sims)
  GxE_mse_vec_pv4 <- numeric(num_sims)

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

    ## 80% of the sample size
    smp_size <- floor(0.8 * (n_e0 + n_e1))

    train_ind <- sample(seq_len(n_e0 + n_e1), size = smp_size)

    train_df <- polygenic_df[train_ind, ]
    test_df <- polygenic_df[-train_ind, ]

    snp_names <- paste0("g", 1:n_snps)

    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    additive_models <- estimate_additive_models(train_df)
    GxE_models <- estimate_GxE_models(train_df)

    additive_ests_pv1 <- get_ests_additive(additive_models, 1)
    GxE_ests_pv1 <- get_ests_GxE(GxE_models, 1)

    additive_ests_pv2 <- get_ests_additive(additive_models, .05)
    GxE_ests_pv2 <- get_ests_GxE(GxE_models, .05)

    additive_ests_pv3 <- get_ests_additive(additive_models, 1e-05)
    GxE_ests_pv3 <- get_ests_GxE(GxE_models, 1e-05)

    additive_ests_pv4 <- get_ests_additive(additive_models, 1e-08)
    GxE_ests_pv4 <- get_ests_GxE(GxE_models, 1e-08)

    additive_preds_pv1 <- test_mat %*% additive_ests_pv1[, 1]
    GxE_preds_pv1 <- c(test_mat_e0 %*% GxE_ests_pv1[, 1], test_mat_e1 %*% GxE_ests_pv1[, 2])

    additive_preds_pv2 <- test_mat %*% additive_ests_pv2[, 1]
    GxE_preds_pv2 <- c(test_mat_e0 %*% GxE_ests_pv2[, 1], test_mat_e1 %*% GxE_ests_pv2[, 2])

    additive_preds_pv3 <- test_mat %*% additive_ests_pv3[, 1]
    GxE_preds_pv3 <- c(test_mat_e0 %*% GxE_ests_pv3[, 1], test_mat_e1 %*% GxE_ests_pv3[, 2])

    additive_preds_pv4 <- test_mat %*% additive_ests_pv4[, 1]
    GxE_preds_pv4 <- c(test_mat_e0 %*% GxE_ests_pv4[, 1], test_mat_e1 %*% GxE_ests_pv4[, 2])

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    #mash_preds <- get_test_preds_mash(Bhat, Shat, test_df)

    # Add metrics to vectors
    #mash_corr_vec[i] <- cor(mash_preds$preds, test_df$y)
    #mash_mse_vec[i] <- mse(
    #  matrix(data = c(mash_preds$est_e0, mash_preds$est_e1), ncol = 2), fx_mat
    #)

    additive_corr_vec_pv1[i] <- protected_cor(additive_preds_pv1, test_df$y)
    GxE_corr_vec_pv1[i] <- protected_cor(GxE_preds_pv1, test_df$y)

    additive_mse_vec_pv1[i] <- mse(additive_ests_pv1, fx_mat)
    GxE_mse_vec_pv1[i] <- mse(GxE_ests_pv1, fx_mat)

    additive_corr_vec_pv2[i] <- protected_cor(additive_preds_pv2, test_df$y)
    GxE_corr_vec_pv2[i] <- protected_cor(GxE_preds_pv2, test_df$y)

    additive_mse_vec_pv2[i] <- mse(additive_ests_pv2, fx_mat)
    GxE_mse_vec_pv2[i] <- mse(GxE_ests_pv2, fx_mat)

    additive_corr_vec_pv3[i] <- protected_cor(additive_preds_pv3, test_df$y)
    GxE_corr_vec_pv3[i] <- protected_cor(GxE_preds_pv3, test_df$y)

    additive_mse_vec_pv3[i] <- mse(additive_ests_pv3, fx_mat)
    GxE_mse_vec_pv3[i] <- mse(GxE_ests_pv3, fx_mat)

    additive_corr_vec_pv4[i] <- protected_cor(additive_preds_pv4, test_df$y)
    GxE_corr_vec_pv4[i] <- protected_cor(GxE_preds_pv4, test_df$y)

    additive_mse_vec_pv4[i] <- mse(additive_ests_pv4, fx_mat)
    GxE_mse_vec_pv4[i] <- mse(GxE_ests_pv4, fx_mat)

  }

  results_list <- list(
   # mash_corr = mean(mash_corr_vec),
  #  mash_mse = mean(mash_mse_vec),

    additive_corr_pv1 = mean(additive_corr_vec_pv1),
    GxE_corr_pv1 = mean(GxE_corr_vec_pv1),

    additive_mse_pv1 = mean(additive_mse_vec_pv1),
    GxE_mse_pv1 = mean(GxE_mse_vec_pv1),

    additive_corr_pv2 = mean(additive_corr_vec_pv2),
    GxE_corr_pv2 = mean(GxE_corr_vec_pv2),

    additive_mse_pv2 = mean(additive_mse_vec_pv2),
    GxE_mse_vec_pv2 = mean(GxE_mse_vec_pv2),

    additive_corr_pv3 = mean(additive_corr_vec_pv3),
    GxE_corr_pv3 = mean(GxE_corr_vec_pv3),

    additive_mse_pv3 = mean(additive_mse_vec_pv3),
    GxE_mse_pv3 = mean(GxE_mse_vec_pv3),

    additive_corr_pv4 = mean(additive_corr_vec_pv4),
    GxE_corr_pv4 = mean(GxE_corr_vec_pv4),

    additive_mse_pv4 = mean(additive_mse_vec_pv4),
    GxE_mse_pv4 = mean(GxE_mse_vec_pv4)
  )

  return(results_list)

}

estimate_GxE_model_fast <- function(n, maf_vec, fx_vec, h2) {

  n_snps <- length(maf_vec)
  var_G <- sum((fx_vec ^ 2) * (2 * (1 - maf_vec) * maf_vec)) - (fx_vec ^ 2) * (2 * (1 - maf_vec) * maf_vec)
  var_E <- get_theo_noise_from_h2(maf_vec, fx_vec, h2) ^ 2
  beta_sd <- sqrt((var_G + var_E) / (2 * (n + 1) * maf_vec * (1 - maf_vec)))
  beta_hat <- rnorm(n_snps, fx_vec, beta_sd)
  shat <- sqrt(((beta_sd ^ 2) / (n - 1)) * rchisq(n_snps, n - 1))
  return(
    list(
      bhat = beta_hat,
      shat = beta_sd,
      pval = pt(-abs(beta_hat / beta_sd), df = n - 2) + (1 - pt(abs(beta_hat / beta_sd), df = n - 2))
    )
  )

}

select_snps_pval_thresh <- function(GxE_models, additive_models, pval = 5e-8) {

  GxE_e0_selection <- GxE_models$pval_e0 < pval
  GxE_e1_selection <- GxE_models$pval_e1 < pval

  additive_selection <- additive_models$pval < pval

  return(
    list(
      e0_selection = GxE_e0_selection,
      e1_selection = GxE_e1_selection,
      additive_selection = additive_selection
    )
  )

}

select_snps_mash <- function(mash_model, lfsr) {

  e0_selection <- mash_model$lfsr_e0 < lfsr
  e1_selection <- mash_model$lfsr_e1 < lfsr

  return(
    list(
      e0_selection = e0_selection,
      e1_selection = e1_selection,
      additive_selection = (e0_selection | e1_selection)
    )
  )

}

select_snps_same_number <- function(GxE_models, additive_models, pval = 5e-8) {

  additive_selection <- additive_models$pval < pval
  num_selected <- sum(additive_selection)

  # Now, want to find the p-value so that I get num_selected for GxE
  sorted_e0 <- sort(GxE_models$pval_e0)
  sorted_e1 <- sort(GxE_models$pval_e1)

  e0_pval <- sorted_e0[num_selected]
  e1_pval <- sorted_e1[num_selected]

  GxE_e0_selection <- GxE_models$pval_e0 <= e0_pval
  GxE_e1_selection <- GxE_models$pval_e1 <= e1_pval

  return(
    list(
      e0_selection = GxE_e0_selection,
      e1_selection = GxE_e1_selection,
      additive_selection = additive_selection
    )
  )

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
simulate_pgs_corr_fast <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  selection_methods = c("additive", "GxE", "mash"),
  beta_methods = c("additive", "GxE", "mash")
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (selection_method in selection_methods) {

    sim_results[[selection_method]] <- list()

    for (model in c("additive", "GxE")) {

      sim_results[[selection_method]][[model]] <- list()

      for (metric in c("mse", "power", "type_1_error", "accuracy", "corr")) {

        sim_results[[selection_method]][[model]][[metric]] <- 0

      }

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    n_e0_train <- ceiling(.8 * n_e0)
    n_e0_test <- n_e0 - n_e0_train

    n_e1_train <- ceiling(.8 * n_e1)
    n_e1_test <- n_e1 - n_e1_train

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic test data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0_test,
      n_e1_test,
      e0_h2,
      e1_h2,
      maf_vec
    )

    test_df <- polygenic_df
    #
    snp_names <- paste0("g", 1:n_snps)
    #
    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()



    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    if ("mash" %in% selection_methods) {

      mash_model <- get_test_preds_mash(Bhat, Shat, test_df)

    }

    for (selection_method in selection_methods) {

      for (beta_method in selection_methods) {

        est_mat_GxE <- NULL

        if (selection_method == "pval") {

          selections <- select_snps_pval_thresh(GxE_models, additive_models)

        } else if (selection_method == "same_number") {

          selections <- select_snps_same_number(GxE_models, additive_models)

        } else if (selection_method == "mash") {

          selections <- select_snps_mash(mash_model)
          est_mat_GxE <- matrix(
            data = c(mash_model$est_e0, mash_model$est_e1),
            ncol = 2
          )

        }

        est_mat_additive <- matrix(
          data = rep(additive_models$est, 2),
          ncol = 2
        )

        sel_mat_additive <- matrix(
          data = rep(selections$additive_selection, 2),
          ncol = 2
        )

        sel_mat_GxE <- matrix(
          data = c(selections$e0_selection, selections$e1_selection),
          ncol = 2
        )

        if (is.null(est_mat_GxE)) {

          est_mat_GxE <- Bhat

        }

        additive_sel_results <- evaluate_selection(
          fx_mat, est_mat_additive, sel_mat_additive
        )

        GxE_sel_results <- evaluate_selection(
          fx_mat, est_mat_GxE, sel_mat_GxE
        )

        additive_selection_ests <- est_mat_additive
        additive_selection_ests[!sel_mat_additive] <- 0
        additive_preds <- test_mat %*% additive_selection_ests[,1]
        additive_sel_results[['corr']] <- protected_cor(additive_preds, test_df$y)

        GxE_selection_ests <- est_mat_GxE
        GxE_selection_ests[!sel_mat_GxE] <- 0
        GxE_preds <- c(
          test_mat_e0 %*% GxE_selection_ests[, 1],
          test_mat_e1 %*% GxE_selection_ests[, 2]
        )

        GxE_sel_results[['corr']] <- protected_cor(GxE_preds, test_df$y)

        comb_sel_results <- list()
        comb_sel_results[["additive"]] <- additive_sel_results
        comb_sel_results[["GxE"]] <- GxE_sel_results

        for (model in c("GxE", "additive")) {

          for (metric in c("mse", "power", "type_1_error", "accuracy", "corr")) {

            sim_results[[selection_method]][[model]][[metric]] <- sim_results[[selection_method]][[model]][[metric]] +
              (1 / num_sims) * comb_sel_results[[model]][[metric]]

          }

        }

      }

    }

  }

  return(sim_results)

}

get_selection_metrics <- function(fx_mat, selection_mat, est_mat) {

  true_signals <- fx_mat != 0
  est_signals <- selection_mat != 0

  power <- sum(true_signals & est_signals) / sum(true_signals)
  type_1_error <- sum(!true_signals & est_signals) / sum(!true_signals)
  accuracy <- sum(est_signals == true_signals) / (2 * nrow(true_signals))

  mse <- mean((est_mat - fx_mat) ^ 2)

  perc_selected <- sum(est_signals) / (2 * nrow(true_signals))

  selected_est_mat <- est_mat
  selected_est_mat[!selection_mat] <- 0
  mse_total <- mean((selected_est_mat - fx_mat) ^ 2)

  perc_correct_sign <- mean(
    sign(est_mat[true_signals & est_signals]) == sign(fx_mat[true_signals & est_signals])
  )

  return(
    list(
      power = power,
      type_1_error = type_1_error,
      accuracy = accuracy,
      mse = mse,
      perc_selected = perc_selected,
      mse_total = mse_total,
      perc_correct_sign = perc_correct_sign
    )
  )

}

get_selection_metrics_bucketed <- function(fx_mat, selection_mat, est_mat) {

  true_bucket_signals <- sum(fx_mat[selection_mat] != 0)
  num_snps_in_bucket <- sum(selection_mat)
  null_snp_pct <- (num_snps_in_bucket - true_bucket_signals) / num_snps_in_bucket

  non_null_fx <- (fx_mat != 0)
  non_null_fx_bucket <- fx_mat[selection_mat & non_null_fx]
  non_null_fx_est <- est_mat[selection_mat & non_null_fx]
  pct_correct_sign <- mean(sign(non_null_fx_bucket) == sign(non_null_fx_est))

  return(
    list(
      null_pct = null_snp_pct,
      correct_sign_pct = pct_correct_sign
    )
  )

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
simulate_pgs_corr_fast_v2 <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  selection_methods = c("additive", "GxE", "mash"),
  pval_thresh_additive = 5e-8,
  pval_thresh_GxE = 5e-8
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (selection_method in selection_methods) {

    sim_results[[selection_method]] <- list()

    for (metric in c("power", "type_1_error", "accuracy", "mse")) {

      sim_results[[selection_method]][[metric]] <- 0

    }

    for (beta_method in selection_methods) {

      sim_results[[selection_method]][[beta_method]] <- list()

      sim_results[[selection_method]][[beta_method]][["corr"]] <- 0

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0,
      n_e1,
      e0_h2,
      e1_h2,
      maf_vec
    )

    ## 80% of the sample size
    smp_size <- floor(0.75 * (n_e0 + n_e1))

    train_ind <- sample(seq_len(n_e0 + n_e1), size = smp_size)

    train_df <- polygenic_df[train_ind, ]
    test_df <- polygenic_df[-train_ind, ]

    snp_names <- paste0("g", 1:n_snps)

    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    # here, I want to estimate the actual GWAS
    additive_models <- estimate_additive_models(train_df)
    GxE_models <- estimate_GxE_models(train_df)

    selected_snps <- list()
    beta_ests <- list()

    selected_snps[["GxE"]] <- matrix(
      data = c(GxE_models$pval_e0 < pval_thresh_GxE, GxE_models$pval_e1 < pval_thresh_GxE),
      ncol = 2
    )

    selected_snps[["additive"]] <- matrix(
      data = rep(additive_models$pval < pval_thresh_additive, 2), ncol = 2
    )

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models$est, 2), ncol = 2
    )

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    beta_ests[['GxE']] <- Bhat

    if ("mash" %in% selection_methods) {

      mash_model <- get_test_preds_mash(
        Bhat, Shat, test_df, amp_range = c(1, 3)
      )
      lfsr_thresh <- 1

      selected_snps[["mash"]] <- matrix(
        data = c(mash_model$lfsr_e0 <= lfsr_thresh, mash_model$lfsr_e1 <= lfsr_thresh),
        ncol = 2
      )

      beta_ests[['mash']] <- matrix(
        data = c(mash_model$est_e0, mash_model$est_e1),
        ncol = 2
      )

    }

    for (selection_method in selection_methods) {

      selection_metrics <- get_selection_metrics(
        fx_mat, selected_snps[[selection_method]], beta_ests[[selection_method]]
      )

      for (metric in c("power", "type_1_error", "accuracy", "mse")) {

        sim_results[[selection_method]][[metric]] <- sim_results[[selection_method]][[metric]] +
          (1 / num_sims) * selection_metrics[[metric]]

      }

      for (beta_method in selection_methods) {

        # get the correlation for each of the prediction methods
        selected_beta_ests <- beta_ests[[beta_method]]
        selected_beta_ests[!selected_snps[[selection_method]]] <- 0

        preds <- c(
          test_mat_e0 %*% selected_beta_ests[, 1],
          test_mat_e1 %*% selected_beta_ests[, 2]
        )

        # NOTE: finish this. Need correct update for simulation object.
        sim_results[[selection_method]][[beta_method]][['corr']] <- sim_results[[selection_method]][[beta_method]][['corr']] +
          (1 / num_sims) * protected_cor(preds, test_df$y)

      }

    }

  }

  return(sim_results)

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
simulate_two_part_mash_sim <- function(
  n_e0_train,
  n_e1_train,
  n_e0_test,
  n_e1_test,
  n_snps_train,
  n_snps_test,
  maf_simulator_train,
  maf_simulator_test,
  Sigma_train,
  Sigma_test,
  pi_train,
  pi_test,
  s_train,
  s_test,
  e0_h2,
  e1_h2,
  num_sims = 10,
  selection_methods = c("additive", "GxE", "mash"),
  pval_thresh_additive = 5e-8,
  pval_thresh_GxE = 5e-8
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (selection_method in selection_methods) {

    sim_results[[selection_method]] <- list()

    for (metric in c("power", "type_1_error", "accuracy", "mse")) {

      sim_results[[selection_method]][[metric]] <- 0

    }

    for (beta_method in selection_methods) {

      sim_results[[selection_method]][[beta_method]] <- list()

      sim_results[[selection_method]][[beta_method]][["corr"]] <- 0

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat_train <- generate_mash_fx(
      n = n_snps_train, p = 2, pi = pi_train, Sigma = Sigma_train, s = s_train
    )

    fx_mat_test <- generate_mash_fx(
      n = n_snps_test, p = 2, pi = pi_test, Sigma = Sigma_test, s = s_test
    )

    # simulate minor allele frequencies all at once for efficiency
    maf_vec_train <- do.call(maf_simulator_train, list(n_snps_train))

    # obtain polygenic data
    polygenic_df_train <- simulate_polygenic_data_2env(
      fx_mat_train, # true effects
      n_e0_train,
      n_e1_train,
      e0_h2,
      e1_h2,
      maf_vec_train
    )

    train_df <- polygenic_df_train

    maf_vec_test <- do.call(maf_simulator_test, list(n_snps_test))

    polygenic_df_test <- simulate_polygenic_data_2env(
      fx_mat_test, # true effects
      n_e0_test,
      n_e1_test,
      e0_h2,
      e1_h2,
      maf_vec_test
    )

    ## 80% of the sample size
    smp_size <- floor(0.75 * (n_e0_test + n_e1_test))

    test_train_ind <- sample(seq_len(n_e0_test + n_e1_test), size = smp_size)

    test_train_df <- polygenic_df_test[test_train_ind, ]
    test_test_df <- polygenic_df_test[-test_train_ind, ]

    snp_names <- paste0("g", 1:n_snps_test)

    test_test_mat <- test_test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_test_mat_e0 <- test_test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_test_mat_e1 <- test_test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    # Estimate models on training and test sets
    GxE_models_train <- estimate_GxE_models(train_df)
    GxE_models_test <- estimate_GxE_models(test_train_df)
    additive_models_test <- estimate_additive_models(test_train_df)

    selected_snps <- list()
    beta_ests <- list()

    selected_snps[["GxE"]] <- matrix(
      data = c(
        GxE_models_test$pval_e0 < pval_thresh_GxE,
        GxE_models_test$pval_e1 < pval_thresh_GxE
      ),
      ncol = 2
    )

    selected_snps[["additive"]] <- matrix(
      data = rep(additive_models_test$pval < pval_thresh_additive, 2), ncol = 2
    )

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models_test$est, 2), ncol = 2
    )

    Bhat_train <- matrix(
      data = c(GxE_models_train$est_e0, GxE_models_train$est_e1), ncol = 2
    )
    Shat_train <- matrix(
      data = c(GxE_models_train$sd_e0, GxE_models_train$sd_e1), ncol = 2
    )

    Bhat_test <- matrix(
      data = c(GxE_models_test$est_e0, GxE_models_test$est_e1), ncol = 2
    )
    Shat_test<- matrix(
      data = c(GxE_models_test$sd_e0, GxE_models_test$sd_e1), ncol = 2
    )

    beta_ests[['GxE']] <- Bhat_test

    if ("mash" %in% selection_methods) {

      mash_train <- fit_mash_gwas_2env(
        Bhat = Bhat_train,
        Shat = Shat_train,
        corr_range = seq(from = -1, to = 1, by = (1/3)),
        amp_range = c(1.5, 2)
      )

      mash_test <- get_test_preds_mash(
        Bhat_test, Shat_test, test_test_df, fitted_g = ashr::get_fitted_g(mash_train)
      )
      lfsr_thresh <- 1

      selected_snps[["mash"]] <- matrix(
        data = c(
          mash_test$lfsr_e0 <= lfsr_thresh, mash_test$lfsr_e1 <= lfsr_thresh
        ),
        ncol = 2
      )

      beta_ests[['mash']] <- matrix(
        data = c(mash_test$est_e0, mash_test$est_e1),
        ncol = 2
      )

    }

    for (selection_method in selection_methods) {

      selection_metrics <- get_selection_metrics(
        fx_mat_test, selected_snps[[selection_method]], beta_ests[[selection_method]]
      )

      for (metric in c("power", "type_1_error", "accuracy", "mse")) {

        sim_results[[selection_method]][[metric]] <- sim_results[[selection_method]][[metric]] +
          (1 / num_sims) * selection_metrics[[metric]]

      }

      for (beta_method in selection_methods) {

        # get the correlation for each of the prediction methods
        selected_beta_ests <- beta_ests[[beta_method]]
        selected_beta_ests[!selected_snps[[selection_method]]] <- 0

        preds <- c(
          test_test_mat_e0 %*% selected_beta_ests[, 1],
          test_test_mat_e1 %*% selected_beta_ests[, 2]
        )

        sim_results[[selection_method]][[beta_method]][['corr']] <- sim_results[[selection_method]][[beta_method]][['corr']] +
          (1 / num_sims) * protected_cor(preds, test_test_df$y)

      }

    }

  }

  return(sim_results)

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
simulate_pgs_corr_one_env <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  selection_methods = c("additive", "GxE", "mash"),
  pval_thresh_additive = 5e-8,
  pval_thresh_GxE = 5e-8
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (selection_method in selection_methods) {

    sim_results[[selection_method]] <- list()

    for (metric in c("power", "type_1_error", "accuracy", "mse")) {

      sim_results[[selection_method]][[metric]] <- 0

    }

    for (beta_method in selection_methods) {

      sim_results[[selection_method]][[beta_method]] <- list()

      sim_results[[selection_method]][[beta_method]][["corr"]] <- 0

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0,
      n_e1,
      e0_h2,
      e1_h2,
      maf_vec
    )

    ## 80% of the sample size
    smp_size <- floor(0.75 * (n_e0 + n_e1))

    train_ind <- sample(seq_len(n_e0 + n_e1), size = smp_size)

    train_df <- polygenic_df[train_ind, ]
    test_df <- polygenic_df[-train_ind, ]

    snp_names <- paste0("g", 1:n_snps)

    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_df_e1 <- test_df %>%
      dplyr::filter(e == "e1")

    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    # here, I want to estimate the actual GWAS
    additive_models <- estimate_additive_models(train_df)
    GxE_models <- estimate_GxE_models(train_df)

    selected_snps <- list()
    beta_ests <- list()

    selected_snps[["GxE"]] <- matrix(
      data = c(GxE_models$pval_e0 < pval_thresh_GxE, GxE_models$pval_e1 < pval_thresh_GxE),
      ncol = 2
    )

    selected_snps[["additive"]] <- matrix(
      data = rep(additive_models$pval < pval_thresh_additive, 2), ncol = 2
    )

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models$est, 2), ncol = 2
    )

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    beta_ests[['GxE']] <- Bhat

    if ("mash" %in% selection_methods) {

      mash_model <- get_test_preds_mash(
        Bhat, Shat, test_df, amp_range = c(1, 3)
      )
      lfsr_thresh <- 1

      selected_snps[["mash"]] <- matrix(
        data = c(mash_model$lfsr_e0 <= lfsr_thresh, mash_model$lfsr_e1 <= lfsr_thresh),
        ncol = 2
      )

      beta_ests[['mash']] <- matrix(
        data = c(mash_model$est_e0, mash_model$est_e1),
        ncol = 2
      )

    }

    for (selection_method in selection_methods) {

      selection_metrics <- get_selection_metrics(
        fx_mat, selected_snps[[selection_method]], beta_ests[[selection_method]]
      )

      for (metric in c("power", "type_1_error", "accuracy", "mse")) {

        sim_results[[selection_method]][[metric]] <- sim_results[[selection_method]][[metric]] +
          (1 / num_sims) * selection_metrics[[metric]]

      }

      for (beta_method in selection_methods) {

        # get the correlation for each of the prediction methods
        selected_beta_ests <- beta_ests[[beta_method]]
        selected_beta_ests[!selected_snps[[selection_method]]] <- 0

        preds <- c(
          test_mat_e1 %*% selected_beta_ests[, 2]
        )

        # NOTE: finish this. Need correct update for simulation object.
        sim_results[[selection_method]][[beta_method]][['corr']] <- sim_results[[selection_method]][[beta_method]][['corr']] +
          (1 / num_sims) * protected_cor(preds, test_df_e1$y)

      }

    }

  }

  return(sim_results)

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
simulate_mash_perf_lfsr <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  lfsr_thresh_vec = seq(.01, 1, .01)
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (lfsr_thresh in lfsr_thresh_vec) {

    sim_results[[glue::glue("thresh_{lfsr_thresh}")]] <- list()

    for (metric in c("power", "type_1_error", "accuracy", "mse", "corr", "perc_selected", "mse_preds")) {

      sim_results[[glue::glue("thresh_{lfsr_thresh}")]][[metric]] <- 0

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    n_e0_train <- ceiling(.75 * n_e0)
    n_e0_test <- n_e0 - n_e0_train

    n_e1_train <- ceiling(.75 * n_e1)
    n_e1_test <- n_e1 - n_e1_train

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic test data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0_test,
      n_e1_test,
      e0_h2,
      e1_h2,
      maf_vec
    )

    test_df <- polygenic_df
    #
    snp_names <- paste0("g", 1:n_snps)
    #
    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    selected_snps <- list()

    GxE_e0 <- estimate_GxE_model_fast(n_e0_train, maf_vec, fx_mat[,1], e0_h2)
    GxE_e1 <- estimate_GxE_model_fast(n_e1_train, maf_vec, fx_mat[,2], e1_h2)

    GxE_models <- list(
      est_e0 = GxE_e0$bhat,
      sd_e0 = GxE_e0$shat,
      pval_e0 = GxE_e0$pval,
      est_e1 = GxE_e1$bhat,
      sd_e1 = GxE_e1$shat,
      pval_e1 = GxE_e1$pval
    )

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    mash_model <- get_test_preds_mash(Bhat, Shat, test_df)

    for (lfsr_thresh in lfsr_thresh_vec) {

      selected_snps <- matrix(
        data = c(mash_model$lfsr_e0 <= lfsr_thresh, mash_model$lfsr_e1 <= lfsr_thresh),
        ncol = 2
      )

      beta_ests <- matrix(
        data = c(mash_model$est_e0, mash_model$est_e1),
        ncol = 2
      )

      selection_metrics <- get_selection_metrics(
        fx_mat, selected_snps, beta_ests
      )

      for (metric in c("power", "type_1_error", "accuracy", "mse", "perc_selected")) {

        sim_results[[glue::glue("thresh_{lfsr_thresh}")]][[metric]] <- sim_results[[glue::glue("thresh_{lfsr_thresh}")]][[metric]] +
          (1 / num_sims) * selection_metrics[[metric]]

      }

      selected_beta_ests <- beta_ests
      selected_beta_ests[!selected_snps] <- 0

      preds <- c(
        test_mat_e0 %*% selected_beta_ests[, 1],
        test_mat_e1 %*% selected_beta_ests[, 2]
      )

      sim_results[[glue::glue("thresh_{lfsr_thresh}")]][['corr']] <- sim_results[[glue::glue("thresh_{lfsr_thresh}")]][['corr']] +
        (1 / num_sims) * protected_cor(preds, test_df$y)

      sim_results[[glue::glue("thresh_{lfsr_thresh}")]][['mse_preds']] <- sim_results[[glue::glue("thresh_{lfsr_thresh}")]][['mse_preds']] +
        (1 / num_sims) * mean((preds - test_df$y) ^ 2)

    }

  }

  return(sim_results)

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
simulate_perf_pval <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  pval_thresh_vec = seq(.01, 1, .01)
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (pval_thresh in pval_thresh_vec) {

    sim_results[[glue::glue("thresh_{pval_thresh}")]] <- list()

    for (model in c("additive", "GxE")) {

      sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]] <- list()

      for (metric in c(
        "mse_preds", "mse_total", "power", "type_1_error", "accuracy",
        "mse", "corr", "perc_selected", "perc_correct_sign")
      ) {

        sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]][[metric]] <- 0

      }

    }

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    n_e0_train <- ceiling(.75 * n_e0)
    n_e0_test <- n_e0 - n_e0_train

    n_e1_train <- ceiling(.75 * n_e1)
    n_e1_test <- n_e1 - n_e1_train

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic test data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0_test,
      n_e1_test,
      e0_h2,
      e1_h2,
      maf_vec
    )

    test_df <- polygenic_df
    #
    snp_names <- paste0("g", 1:n_snps)
    #
    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    GxE_e0 <- estimate_GxE_model_fast(n_e0_train, maf_vec, fx_mat[,1], e0_h2)
    GxE_e1 <- estimate_GxE_model_fast(n_e1_train, maf_vec, fx_mat[,2], e1_h2)

    GxE_models <- list(
      est_e0 = GxE_e0$bhat,
      sd_e0 = GxE_e0$shat,
      pval_e0 = GxE_e0$pval,
      est_e1 = GxE_e1$bhat,
      sd_e1 = GxE_e1$shat,
      pval_e1 = GxE_e1$pval
    )

    additive_models_bhat <- (
      (1 / GxE_e0$shat ^ 2) * GxE_e0$bhat + (1 / GxE_e1$shat ^ 2) * GxE_e1$bhat
    ) / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2))

    additive_models_shat <- sqrt(1 / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2)))

    additive_models_pval <- pt(
      -abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2
    ) +
      (1 - pt(abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2))

    additive_models <- list(
      est = additive_models_bhat,
      pval = additive_models_pval
    )

    beta_ests <- list()

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models$est, 2), ncol = 2
    )

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    beta_ests[['GxE']] <- Bhat

    for (pval_thresh in pval_thresh_vec) {

      selected_snps <- list()

      selected_snps[["GxE"]] <- matrix(
        data = c(GxE_models$pval_e0 < pval_thresh, GxE_models$pval_e1 < pval_thresh),
        ncol = 2
      )

      selected_snps[["additive"]] <- matrix(
        data = rep(additive_models$pval < pval_thresh, 2), ncol = 2
      )

      # Here, get selection metrics for the additive and GxE models

      selection_metrics <- list()

      selection_metrics[["additive"]] <- get_selection_metrics(
        fx_mat, selected_snps[["additive"]], beta_ests[["additive"]]
      )

      selection_metrics[["GxE"]] <- get_selection_metrics(
        fx_mat, selected_snps[["GxE"]], beta_ests[["GxE"]]
      )

      for (metric in c("mse_total", "power", "type_1_error", "accuracy", "mse",
                       "perc_selected", "perc_correct_sign")) {

        for (model in c("additive", "GxE")) {

          sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]][[metric]] <- sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]][[metric]] +
            (1 / num_sims) * selection_metrics[[model]][[metric]]

        }

      }

      # Now, want an additive and GxE corr

      selected_beta_ests_additive <- beta_ests[['additive']]
      selected_beta_ests_additive[!selected_snps[["additive"]]] <- 0

      additive_preds <- c(
        test_mat_e0 %*% selected_beta_ests_additive[, 1],
        test_mat_e1 %*% selected_beta_ests_additive[, 2]
      )

      sim_results[[glue::glue("thresh_{pval_thresh}")]][["additive"]][['corr']] <- sim_results[[glue::glue("thresh_{pval_thresh}")]][["additive"]][['corr']] +
        (1 / num_sims) * protected_cor(additive_preds, test_df$y)

      sim_results[[glue::glue("thresh_{pval_thresh}")]][["additive"]][['mse_preds']] <- sim_results[[glue::glue("thresh_{pval_thresh}")]][["additive"]][['mse_preds']] +
        (1 / num_sims) * mean((additive_preds - test_df$y) ^ 2)

      selected_beta_ests_GxE <- beta_ests[['GxE']]
      selected_beta_ests_GxE[!selected_snps[["GxE"]]] <- 0

      GxE_preds <- c(
        test_mat_e0 %*% selected_beta_ests_GxE[, 1],
        test_mat_e1 %*% selected_beta_ests_GxE[, 2]
      )

      sim_results[[glue::glue("thresh_{pval_thresh}")]][["GxE"]][['corr']] <- sim_results[[glue::glue("thresh_{pval_thresh}")]][["GxE"]][['corr']] +
        (1 / num_sims) * protected_cor(GxE_preds, test_df$y)

      sim_results[[glue::glue("thresh_{pval_thresh}")]][["GxE"]][['mse_preds']] <- sim_results[[glue::glue("thresh_{pval_thresh}")]][["GxE"]][['mse_preds']] +
        (1 / num_sims) * mean((GxE_preds - test_df$y) ^ 2)

    }

  }

  return(sim_results)

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
simulate_pval_buckets <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  pval_thresh_vec = seq(.01, 1, .01)
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (pval_thresh in pval_thresh_vec) {

    sim_results[[glue::glue("thresh_{pval_thresh}")]] <- list()

    for (model in c("additive", "GxE")) {

      sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]] <- list()

      for (metric in c(
        "null_pct", "correct_sign_pct")
      ) {

        sim_results[[glue::glue("thresh_{pval_thresh}")]][[model]][[metric]] <- 0

      }

    }

  }

  pval_thresh_vec <- c(0, pval_thresh_vec)

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    n_e0_train <- ceiling(.75 * n_e0)
    n_e0_test <- n_e0 - n_e0_train

    n_e1_train <- ceiling(.75 * n_e1)
    n_e1_test <- n_e1 - n_e1_train

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    GxE_e0 <- estimate_GxE_model_fast(n_e0_train, maf_vec, fx_mat[,1], e0_h2)
    GxE_e1 <- estimate_GxE_model_fast(n_e1_train, maf_vec, fx_mat[,2], e1_h2)

    GxE_models <- list(
      est_e0 = GxE_e0$bhat,
      sd_e0 = GxE_e0$shat,
      pval_e0 = GxE_e0$pval,
      est_e1 = GxE_e1$bhat,
      sd_e1 = GxE_e1$shat,
      pval_e1 = GxE_e1$pval
    )

    additive_models_bhat <- (
      (1 / GxE_e0$shat ^ 2) * GxE_e0$bhat + (1 / GxE_e1$shat ^ 2) * GxE_e1$bhat
    ) / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2))

    additive_models_shat <- sqrt(1 / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2)))

    additive_models_pval <- pt(
      -abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2
    ) +
      (1 - pt(abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2))

    additive_models <- list(
      est = additive_models_bhat,
      pval = additive_models_pval
    )

    beta_ests <- list()

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models$est, 2), ncol = 2
    )

    Bhat <- matrix(
      data = c(GxE_models$est_e0, GxE_models$est_e1), ncol = 2
    )
    Shat <- matrix(
      data = c(GxE_models$sd_e0, GxE_models$sd_e1), ncol = 2
    )

    beta_ests[['GxE']] <- Bhat

    for (i in seq(2, length(pval_thresh_vec))) {

      # here, calculate each of the things that I need above

      selected_snps <- list()

      selected_snps[["GxE"]] <- matrix(
        data = c(
          GxE_models$pval_e0 < pval_thresh_vec[i] &
            GxE_models$pval_e0 >= pval_thresh_vec[i - 1],
          GxE_models$pval_e1 < pval_thresh_vec[i] &
            GxE_models$pval_e1 >= pval_thresh_vec[i - 1]
        ),
        ncol = 2
      )

      selected_snps[["additive"]] <- matrix(
        data = rep(
          additive_models$pval < pval_thresh_vec[i] &
            additive_models$pval >= pval_thresh_vec[i - 1],
          2
        ),
        ncol = 2
      )

      # Here, get selection metrics for the additive and GxE models

      selection_metrics <- list()

      # next, I have to implement the bucketed version of selection metrics
      selection_metrics[["additive"]] <- get_selection_metrics_bucketed(
        fx_mat, selected_snps[["additive"]], beta_ests[["additive"]]
      )

      selection_metrics[["GxE"]] <- get_selection_metrics_bucketed(
        fx_mat, selected_snps[["GxE"]], beta_ests[["GxE"]]
      )

      for (metric in c("null_pct", "correct_sign_pct")) {

        for (model in c("additive", "GxE")) {

          sim_results[[glue::glue("thresh_{pval_thresh_vec[i]}")]][[model]][[metric]] <- sim_results[[glue::glue("thresh_{pval_thresh_vec[i]}")]][[model]][[metric]] +
            (1 / num_sims) * selection_metrics[[model]][[metric]]

        }

      }

    }

  }

  return(sim_results)

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
simulate_perfect_annealing <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  num_sims = 10,
  s = 0,
  perc_included_vec = seq(.01, 1, length.out = 250)
) {

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  sim_results <- list()

  for (perc_included in perc_included_vec) {

    sim_results[[glue::glue("perc_included_{perc_included}")]] <- list()
    sim_results[[glue::glue("perc_included_{perc_included}")]][["corr"]] <- 0
    sim_results[[glue::glue("perc_included_{perc_included}")]][["mse_preds"]] <- 0

  }

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    n_e0_train <- ceiling(.75 * n_e0)
    n_e0_test <- n_e0 - n_e0_train

    n_e1_train <- ceiling(.75 * n_e1)
    n_e1_test <- n_e1 - n_e1_train

    # simulate minor allele frequencies all at once for efficiency
    maf_vec <- do.call(maf_simulator, list(n_snps))

    # obtain polygenic test data
    polygenic_df <- simulate_polygenic_data_2env(
      fx_mat, # true effects
      n_e0_test,
      n_e1_test,
      e0_h2,
      e1_h2,
      maf_vec
    )

    test_df <- polygenic_df
    #
    snp_names <- paste0("g", 1:n_snps)
    #
    test_mat <- test_df %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e0 <- test_df %>%
      dplyr::filter(e == "e0") %>%
      dplyr::select(snp_names) %>%
      as.matrix()
    #
    test_mat_e1 <- test_df %>%
      dplyr::filter(e == "e1") %>%
      dplyr::select(snp_names) %>%
      as.matrix()

    GxE_e0 <- estimate_GxE_model_fast(n_e0_train, maf_vec, fx_mat[,1], e0_h2)
    GxE_e1 <- estimate_GxE_model_fast(n_e1_train, maf_vec, fx_mat[,2], e1_h2)

    GxE_models <- list(
      est_e0 = GxE_e0$bhat,
      sd_e0 = GxE_e0$shat,
      pval_e0 = GxE_e0$pval,
      est_e1 = GxE_e1$bhat,
      sd_e1 = GxE_e1$shat,
      pval_e1 = GxE_e1$pval
    )

    additive_models_bhat <- (
      (1 / GxE_e0$shat ^ 2) * GxE_e0$bhat + (1 / GxE_e1$shat ^ 2) * GxE_e1$bhat
    ) / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2))

    additive_models_shat <- sqrt(1 / ((1 / GxE_e0$shat ^ 2) + (1 / GxE_e1$shat ^ 2)))

    additive_models_pval <- pt(
      -abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2
    ) +
      (1 - pt(abs(additive_models_bhat / additive_models_shat), df = n_e0_train + n_e1_train - 2))

    additive_models <- list(
      est = additive_models_bhat,
      pval = additive_models_pval
    )

    beta_ests <- list()

    beta_ests[['additive']] <- matrix(
      data = rep(additive_models$est, 2), ncol = 2
    )

    additive_pvals_mat <- matrix(
      data = rep(additive_models$pval, 2), ncol = 2
    )

    true_signals <- fx_mat != 0

    sorted_pvals_true <- sort(additive_models_pval[true_signals])
    sorted_pvals_false <- sort(additive_models_pval[!true_signals])

    for (perc_included in perc_included_vec) {

      selected_snps <- list()

      if (perc_included <= (1 - s)) {

        pval_thresh <- sorted_pvals_true[
          min(floor(perc_included * n_snps), floor(n_snps * (1 - s)))
        ]

        additive_pvals_tresh_mat <- (additive_pvals_mat <= pval_thresh)
        selected_snps <- true_signals & additive_pvals_tresh_mat

      } else {

        perc_false_signals <- perc_included - (1 - s)
        pval_thresh <- sorted_pvals_false[
          min(floor(perc_false_signals * n_snps), floor(n_snps * s))
        ]
        additive_pvals_tresh_mat <- (additive_pvals_mat <= pval_thresh)
        selected_snps <- true_signals | additive_pvals_tresh_mat

      }

      # Now, want an additive and GxE corr

      selected_beta_ests_additive <- beta_ests[['additive']]
      selected_beta_ests_additive[!selected_snps] <- 0

      additive_preds <- c(
        test_mat_e0 %*% selected_beta_ests_additive[, 1],
        test_mat_e1 %*% selected_beta_ests_additive[, 2]
      )

      sim_results[[glue::glue("perc_included_{perc_included}")]][['corr']] <- sim_results[[glue::glue("perc_included_{perc_included}")]][['corr']] +
        (1 / num_sims) * protected_cor(additive_preds, test_df$y)

      sim_results[[glue::glue("perc_included_{perc_included}")]][['mse_preds']] <- sim_results[[glue::glue("perc_included_{perc_included}")]][['mse_preds']] +
        (1 / num_sims) * mean((additive_preds - test_df$y) ^ 2)

    }

  }

  return(sim_results)

}

#' Simulate ascertainment from a mixture model
#'
#' @param n_e0,n_e1 number of individuals in each environment
#' @param n_snps number of SNPs
#' @param maf_simulator function that takes number of snps as an argument and
#' returns minor allele frequencies for each snp
#' @param sigma_e0,sigma_e1 additive noise standard deviation in each environment
#' @param Sigma list of covariance matrices from which true effects are drawn
#' @param pi vector of probabilities corresponding to Sigma
#' @param num_sims number of simulations to run
#' @param pval_thresh threshold for p-value
#' @param s (optional) level of sparsity of true effects. Defaults to 0
#'
#' @return
#' @export
#'
#' @examples
simulate_ascertainment <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  pval_thresh,
  num_sims = 10,
  s = 0
) {

  additive_power_vec <- numeric(num_sims)
  GxE_power_vec <- numeric(num_sims)

  additive_tp_mse_vec <- numeric(num_sims)
  GxE_tp_mse_vec <- numeric(num_sims)

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

    train_df <- polygenic_df

    snp_names <- paste0("g", 1:n_snps)

    additive_models <- estimate_additive_models(train_df)
    GxE_models <- estimate_GxE_models(train_df)

    additive_ests <- get_ests_additive(additive_models, pval_thresh)
    GxE_ests <- get_ests_GxE(GxE_models, pval_thresh)

    additive_asc <- get_asc(fx_mat, additive_ests)
    GxE_asc <- get_asc(fx_mat, GxE_ests)

    additive_power_vec[i] <- additive_asc$power
    additive_tp_mse_vec[i] <- additive_asc$tp_mse

    GxE_power_vec[i] <- GxE_asc$power
    GxE_tp_mse_vec[i] <- GxE_asc$tp_mse

  }

  results_list <- list(

    additive_power = mean(additive_power_vec),
    additive_tp_mse = mean(additive_tp_mse_vec),

    GxE_power = mean(GxE_power_vec),
    GxE_tp_mse = mean(GxE_tp_mse_vec)

  )

  return(results_list)

}


#' Simulate ascertainment from a mixture model
#'
#' @param n_e0,n_e1 number of individuals in each environment
#' @param n_snps number of SNPs
#' @param maf_simulator function that takes number of snps as an argument and
#' returns minor allele frequencies for each snp
#' @param sigma_e0,sigma_e1 additive noise standard deviation in each environment
#' @param Sigma list of covariance matrices from which true effects are drawn
#' @param pi vector of probabilities corresponding to Sigma
#' @param num_sims number of simulations to run
#' @param pval_thresh threshold for p-value
#' @param s (optional) level of sparsity of true effects. Defaults to 0
#'
#' @return
#' @export
#'
#' @examples
simulate_ascertainment_fast <- function(
  n_e0,
  n_e1,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  pval_thresh,
  num_sims = 10,
  s = 0
) {

  additive_power_vec <- numeric(num_sims)
  GxE_power_vec <- numeric(num_sims)

  additive_tp_mse_vec <- numeric(num_sims)
  GxE_tp_mse_vec <- numeric(num_sims)

  pb <- txtProgressBar(min = 0, max = num_sims, initial = 0, style = 3)

  for (i in 1:num_sims) {

    setTxtProgressBar(pb, i)

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

    maf_vec <- maf_simulator(n_snps)

    GxE_e0 <- estimate_GxE_model_fast(n_e0, maf_vec, fx_mat[,1], e0_h2)
    GxE_e1 <- estimate_GxE_model_fast(n_e1, maf_vec, fx_mat[,2], e1_h2)

    additive_models_bhat <- (
      (1 / GxE_e0$shat) * GxE_e0$bhat + (1 / GxE_e1$shat) * GxE_e1$bhat
    ) / ((1 / GxE_e0$shat) + (1 / GxE_e1$shat))

    additive_models_shat <- 1 / ((1 / GxE_e0$shat) + (1 / GxE_e1$shat))

    additive_models_pval <- pt(-abs(additive_models_bhat / additive_models_shat), df = n_e0 + n_e1 - 2) +
      (1 - pt(additive_models_bhat / additive_models_shat, df = n_e0 + n_e1 - 2))

    additive_models <- list(
      est = additive_models_bhat,
      pval = additive_models_pval
    )

    GxE_models <- list(
      est_e0 = GxE_e0$bhat,
      sd_e0 = GxE_e0$shat,
      pval_e0 = GxE_e0$pval,
      est_e1 = GxE_e1$bhat,
      sd_e1 = GxE_e1$shat,
      pval_e1 = GxE_e1$pval
    )

    additive_ests <- get_ests_additive(additive_models, pval_thresh)
    GxE_ests <- get_ests_GxE(GxE_models, pval_thresh)

    additive_asc <- get_asc(fx_mat, additive_ests)
    GxE_asc <- get_asc(fx_mat, GxE_ests)

    additive_power_vec[i] <- additive_asc$power
    additive_tp_mse_vec[i] <- additive_asc$tp_mse

    GxE_power_vec[i] <- GxE_asc$power
    GxE_tp_mse_vec[i] <- GxE_asc$tp_mse

  }

  results_list <- list(

    additive_power = mean(additive_power_vec),
    additive_tp_mse = mean(additive_tp_mse_vec),

    GxE_power = mean(GxE_power_vec),
    GxE_tp_mse = mean(GxE_tp_mse_vec)

  )

  return(results_list)

}
