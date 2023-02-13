#' Get Theoretical Observation Noise from GWAS Results
#'
#' This code assumes that y = b * g + e where e ~ N(0, sigma ^ 2), y is a vector
#' of continuous outcomes, and g is a vector of genotypes. Based on an expected
#' standard error for an estimate of b, this function calculates our expected
#' sigma.
#'
#' @param se desired standard error of GWAS effect size.
#' @param n number of individuals in GWAS
#' @param maf minor allele frequency
#'
#' @return expected theoretical observation noise sigma.
#'
get_theo_obs_noise <- function(se, n, maf) {

  sqrt((se ^ 2) * (n + 1) * 2 * maf * (1 - maf))

}

#' Get MSE of GxE model
#'
#' This code assumes we're estimating the model y = b * g + e,
#' where e ~ N(0, sigma ^ 2), separately in two environments.
#'
#' @param sigma_e0,sigma_e1 observation noise in each environment
#' @param g_e0,g_e1 vector of genotypes in each environment
#'
#' @return list with bias, variance, and mse
#'
get_GxE_2env_mse <- function(
  sigma_e0,
  sigma_e1,
  g_e0,
  g_e1
) {

  ssg_e0 <- sum((g_e0 - mean(g_e0)) ^ 2)
  var_e0 <- (sigma_e0 ^ 2) / ssg_e0

  ssg_e1 <- sum((g_e1 - mean(g_e1)) ^ 2)
  var_e1 <- (sigma_e1 ^ 2) / ssg_e1

  return(
    list(
      mse_e0 = var_e0,
      bias_e0 = 0,
      var_e0 = var_e0,
      mse_e1 = var_e1,
      bias_e1 = 0,
      var_e1 = var_e1,
      mse = .5 * var_e0 + .5 * var_e1
    )
  )

}

#' Get MSE For Two Environment Additive Model
#'
#' This function assumes that we're estimating the model y = b * g + e jointly
#' for two environments.
#'
#' @param fx_e0,fx_e1 effect size in each environment
#' @param sigma_e0,sigma_e1 additive noise in each environment
#' @param g_e0,g_e1 vector of genotypes in each environment
#' @param int_e0,int_e1 intercept in each environment. Defaults to 0
#'
#' @return list with bias, variance, and mse.
#'
get_additive_2env_mse <- function(
  fx_e0,
  fx_e1,
  sigma_e0,
  sigma_e1,
  g_e0,
  g_e1,
  int_e0,
  int_e1
) {

  g <- c(g_e0, g_e1)
  g_mean <- mean(g)
  ssg_e0 <- sum((g_e0 - g_mean) ^ 2)
  ssg_e1 <- sum((g_e1 - g_mean) ^ 2)
  ssg <- sum((g - g_mean) ^ 2)

  variance <- ((sigma_e0 ^ 2) * ssg_e0 + (sigma_e1 ^ 2) * ssg_e1) / (ssg ^ 2)

  ssg_slope_e0 <- sum((g_e0 - g_mean) * g_e0)
  ssg_slope_e1 <- sum((g_e1 - g_mean) * g_e1)

  reg_est <- (
    int_e0 * ssg_e0 + int_e1 * ssg_e1 + fx_e0 * ssg_slope_e0 + fx_e1 * ssg_slope_e1
    ) / ssg

  bias_e0 <- reg_est - fx_e0
  bias_e1 <- reg_est - fx_e1

  mse_e0 <- bias_e0 ^ 2 + variance
  mse_e1 <- bias_e1 ^ 2 + variance

  return(
    list(
      mse_e0 = mse_e0,
      bias_e0 = bias_e0,
      var_e0 = variance,
      mse_e1 = mse_e1,
      bias_e1 = bias_e1,
      var_e1 = variance,
      mse = .5 * mse_e0 + .5 * mse_e1
    )
  )

}

#' Simulate MSE for GxE and Additive Model on A Grid With Effect Size Differences
#' And Equal Environmental Variance
#'
#' This code allows the user to specify a grid of standard errors for the
#' GWAS coefficients (where noise is assumed equivalent between environments)
#' and a grid of differences in effect size between environments
#'
#' @param n_e0,n_e1 number of individuals in each environment
#' @param se_grid vector of standard errors
#' @param fx_diff_grid vector of effect size differences between environments
#' @param maf_e0,maf_e1 (optional) minor allele frequencies in each environment.
#' Defaults to .4
#' @param int_e0,int_e1 (optional) intercepts in each environment. Defaults to
#' 0.
#' @param sims_per_grid_pt (optional) number of simulations for each combination
#' of value in \code{se_grid} and \code{fx_norm_diff_grid}. Defaults to 10
#'
#' @return dataframe with various metrics for the GxE and additive models.
#' @export
#'
#' @importFrom magrittr %>%
#'
sim_fx_diff_noise_tradeoff_2env <- function(
  n_e0,
  n_e1,
  se_grid,
  fx_diff_grid,
  maf_e0 = .4,
  maf_e1 = .4,
  int_e0 = 0,
  int_e1 = 0,
  sims_per_grid_pt = 10
) {

  sim_df <- expand.grid(fx_diff_grid, se_grid)
  colnames(sim_df) <- c("fx_diff", "se")

  mse_GxE_vec <- numeric(nrow(sim_df))
  mse_additive_vec <- numeric(nrow(sim_df))

  mse_e0_GxE_vec <- numeric(nrow(sim_df))
  var_e0_GxE_vec <- numeric(nrow(sim_df))
  bias_e0_GxE_vec <- numeric(nrow(sim_df))

  mse_e1_GxE_vec <- numeric(nrow(sim_df))
  var_e1_GxE_vec <- numeric(nrow(sim_df))
  bias_e1_GxE_vec <- numeric(nrow(sim_df))

  mse_e0_additive_vec <- numeric(nrow(sim_df))
  var_e0_additive_vec <- numeric(nrow(sim_df))
  bias_e0_additive_vec <- numeric(nrow(sim_df))

  mse_e1_additive_vec <- numeric(nrow(sim_df))
  var_e1_additive_vec <- numeric(nrow(sim_df))
  bias_e1_additive_vec <- numeric(nrow(sim_df))

  cat("Running simulations...\n")
  pb <- txtProgressBar(min = 0, max = nrow(sim_df), initial = 0, style = 3)

  for (i in 1:nrow(sim_df)) {

    setTxtProgressBar(pb, i)

    mse_sep_iter_vec <- numeric(sims_per_grid_pt)
    mse_comb_iter_vec <- numeric(sims_per_grid_pt)

    sigma_e0 <- get_theo_obs_noise(sim_df$se[i], n_e0, maf_e0)
    sigma_e1 <- get_theo_obs_noise(sim_df$se[i], n_e1, maf_e1)

    for(j in 1:sims_per_grid_pt) {

      g_e0 <- rbinom(n = n_e0, size = 2, prob = maf_e0)
      g_e1 <- rbinom(n = n_e1, size = 2, prob = maf_e1)

      res_GxE <- get_GxE_2env_mse(sigma_e0, sigma_e1, g_e0, g_e1)
      res_additive <- get_additive_2env_mse(
        fx_e0 = 0,
        fx_e1 = sim_df$fx_diff[i],
        int_e0 = int_e0,
        int_e1 = int_e1,
        sigma_e0 = sigma_e0,
        sigma_e1 = sigma_e1,
        g_e0 = g_e0,
        g_e1 = g_e1
      )

      # update vectors
      mse_GxE_vec[i] <- mse_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse
      mse_additive_vec[i] <- mse_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse

      mse_e0_GxE_vec[i] <- mse_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse_e0
      var_e0_GxE_vec[i] <- var_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$var_e0
      bias_e0_GxE_vec[i] <- bias_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$bias_e0

      mse_e1_GxE_vec[i] <- mse_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse_e1
      var_e1_GxE_vec[i] <- var_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$var_e1
      bias_e1_GxE_vec[i] <- bias_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$bias_e1

      mse_e0_additive_vec[i] <- mse_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse_e0
      var_e0_additive_vec[i] <- var_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$var_e0
      bias_e0_additive_vec[i] <- bias_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$bias_e0

      mse_e1_additive_vec[i] <- mse_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse_e1
      var_e1_additive_vec[i] <- var_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$var_e1
      bias_e1_additive_vec[i] <- bias_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$bias_e1

    }

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = mse_GxE_vec,
      mse_additive = mse_additive_vec,

      mse_e0_GxE = mse_e0_GxE_vec,
      var_e0_GxE = var_e0_GxE_vec,
      bias_e0_GxE = bias_e0_GxE_vec,

      mse_e1_GxE = mse_e1_GxE_vec,
      var_e1_GxE = var_e1_GxE_vec,
      bias_e1_GxE = bias_e1_GxE_vec,

      mse_e0_additive = mse_e0_additive_vec,
      var_e0_additive = var_e0_additive_vec,
      bias_e0_additive = bias_e0_additive_vec,

      mse_e1_additive = mse_e1_additive_vec,
      var_e1_additive = var_e1_additive_vec,
      bias_e1_additive = bias_e1_additive_vec
    )

  return(sim_df)

}


#' Simulate MSE for GxE and Additive Model on A Grid With Difference Effect Sizes
#' And Environmental Variance Ratios
#'
#' This code allows the user to specify a grid of ratios of environmental
#' variance between environments
#' and a grid of differences in effect size between environments
#' (normalized by the standard error).
#'
#' @param n_e0,n_e1 number of individuals in each environment
#' @param se_ratio_grid vector of standard error ratios. The standard error
#' in environment 0 will be taken as \code{ref_se} and the error in environment
#' 1 will be set based on this ratio.
#' @param fx_diff_norm_grid vector of effect size differences between environments
#' normalized by trait variance.
#' @param ref_se standard error for environment 0. Used to determine standard
#' error for environment 1.
#' @param maf_e0,maf_e1 (optional) minor allele frequencies in each environment.
#' Defaults to .4
#' @param int_e0,int_e1 (optional) intercepts in each environment. Defaults to
#' 0.
#' @param sims_per_grid_pt (optional) number of simulations for each combination
#' of value in \code{se_grid} and \code{fx_norm_diff_grid}. Defaults to 10
#'
#' @return dataframe with various metrics for the GxE and additive models.
#' @export
#'
#' @importFrom magrittr %>%
#'
sim_diff_variance_2env <- function(
  n_e0,
  n_e1,
  ref_se,
  se_ratio_grid,
  fx_diff_norm_grid,
  maf_e0 = .4,
  maf_e1 = .4,
  int_e0 = 0,
  int_e1 = 0,
  sims_per_grid_pt = 10
) {

  sim_df <- expand.grid(fx_diff_norm_grid, se_ratio_grid)
  colnames(sim_df) <- c("fx_diff_norm", "se_ratio")

  mse_GxE_vec <- numeric(nrow(sim_df))
  mse_additive_vec <- numeric(nrow(sim_df))

  mse_e0_GxE_vec <- numeric(nrow(sim_df))
  var_e0_GxE_vec <- numeric(nrow(sim_df))
  bias_e0_GxE_vec <- numeric(nrow(sim_df))

  mse_e1_GxE_vec <- numeric(nrow(sim_df))
  var_e1_GxE_vec <- numeric(nrow(sim_df))
  bias_e1_GxE_vec <- numeric(nrow(sim_df))

  mse_e0_additive_vec <- numeric(nrow(sim_df))
  var_e0_additive_vec <- numeric(nrow(sim_df))
  bias_e0_additive_vec <- numeric(nrow(sim_df))

  mse_e1_additive_vec <- numeric(nrow(sim_df))
  var_e1_additive_vec <- numeric(nrow(sim_df))
  bias_e1_additive_vec <- numeric(nrow(sim_df))

  cat("Running simulations...\n")
  pb <- txtProgressBar(min = 0, max = nrow(sim_df), initial = 0, style = 3)

  sigma_e0 <- get_theo_obs_noise(ref_se, n_e0, maf_e0)

  for (i in 1:nrow(sim_df)) {

    setTxtProgressBar(pb, i)

    mse_sep_iter_vec <- numeric(sims_per_grid_pt)
    mse_comb_iter_vec <- numeric(sims_per_grid_pt)

    sigma_e1 <- sigma_e0 * sim_df$se_ratio[i]

    implied_fx_e1 <- sim_df$fx_diff_norm[i] * sigma_e0

    for(j in 1:sims_per_grid_pt) {

      g_e0 <- rbinom(n = n_e0, size = 2, prob = maf_e0)
      g_e1 <- rbinom(n = n_e1, size = 2, prob = maf_e1)

      res_GxE <- get_GxE_2env_mse(sigma_e0, sigma_e1, g_e0, g_e1)
      res_additive <- get_additive_2env_mse(
        fx_e0 = 0,
        fx_e1 = implied_fx_e1,
        int_e0 = int_e0,
        int_e1 = int_e1,
        sigma_e0 = sigma_e0,
        sigma_e1 = sigma_e1,
        g_e0 = g_e0,
        g_e1 = g_e1
      )

      # update vectors
      mse_GxE_vec[i] <- mse_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse
      mse_additive_vec[i] <- mse_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse

      mse_e0_GxE_vec[i] <- mse_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse_e0
      var_e0_GxE_vec[i] <- var_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$var_e0
      bias_e0_GxE_vec[i] <- bias_e0_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$bias_e0

      mse_e1_GxE_vec[i] <- mse_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$mse_e1
      var_e1_GxE_vec[i] <- var_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$var_e1
      bias_e1_GxE_vec[i] <- bias_e1_GxE_vec[i] + (1 / sims_per_grid_pt) * res_GxE$bias_e1

      mse_e0_additive_vec[i] <- mse_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse_e0
      var_e0_additive_vec[i] <- var_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$var_e0
      bias_e0_additive_vec[i] <- bias_e0_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$bias_e0

      mse_e1_additive_vec[i] <- mse_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$mse_e1
      var_e1_additive_vec[i] <- var_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$var_e1
      bias_e1_additive_vec[i] <- bias_e1_additive_vec[i] + (1 / sims_per_grid_pt) * res_additive$bias_e1

    }

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = mse_GxE_vec,
      mse_additive = mse_additive_vec,

      mse_e0_GxE = mse_e0_GxE_vec,
      var_e0_GxE = var_e0_GxE_vec,
      bias_e0_GxE = bias_e0_GxE_vec,

      mse_e1_GxE = mse_e1_GxE_vec,
      var_e1_GxE = var_e1_GxE_vec,
      bias_e1_GxE = bias_e1_GxE_vec,

      mse_e0_additive = mse_e0_additive_vec,
      var_e0_additive = var_e0_additive_vec,
      bias_e0_additive = bias_e0_additive_vec,

      mse_e1_additive = mse_e1_additive_vec,
      var_e1_additive = var_e1_additive_vec,
      bias_e1_additive = bias_e1_additive_vec
    )

  return(sim_df)

}

#' Simulate MSE for GxE and Additive Model on A Grid With Difference Effect Sizes
#' And Environmental Variance Ratios
#'
#' This code allows the user to specify a grid of ratios of environmental
#' variance between environments
#' and a grid of differences in effect size between environments
#' (normalized by the standard error).
#'
#' @param sigma_ratio ratio of standard errors between environments
#' @param fx_diff_grid grid of differences in effect size
#' @param se_grid grid of standard errors of regression coefficients
#'
#' @return dataframe with various metrics for the GxE and additive models.
#' @export
#'
#' @importFrom magrittr %>%
#'
sim_var_ratio_2env <- function(
  sigma_ratio,
  fx_diff_grid,
  se_grid
) {

  sim_df <- expand.grid(fx_diff_grid, se_grid)
  colnames(sim_df) <- c("fx_diff", "se")

  mse_GxE_vec <- numeric(nrow(sim_df))
  mse_additive_vec <- numeric(nrow(sim_df))

  mse_e0_GxE_vec <- numeric(nrow(sim_df))
  var_e0_GxE_vec <- numeric(nrow(sim_df))
  bias_e0_GxE_vec <- numeric(nrow(sim_df))

  mse_e1_GxE_vec <- numeric(nrow(sim_df))
  var_e1_GxE_vec <- numeric(nrow(sim_df))
  bias_e1_GxE_vec <- numeric(nrow(sim_df))

  mse_e0_additive_vec <- numeric(nrow(sim_df))
  var_e0_additive_vec <- numeric(nrow(sim_df))
  bias_e0_additive_vec <- numeric(nrow(sim_df))

  mse_e1_additive_vec <- numeric(nrow(sim_df))
  var_e1_additive_vec <- numeric(nrow(sim_df))
  bias_e1_additive_vec <- numeric(nrow(sim_df))

  avg_se_vec <- numeric(nrow(sim_df))

  for (i in 1:nrow(sim_df)) {

    fx_e0 <- 0
    fx_e1 <- sim_df$fx_diff[i]

    se_0 <- sim_df$se[i]
    se_1 <- se_0 * sigma_ratio

    se_additive <- sqrt(.25 * (se_0 ^ 2) + .25 * (se_1 ^ 2))

    avg_se_vec[i] <- se_additive

    exp_additive <- mean(c(fx_e0, fx_e1))

    # update vectors
    mse_GxE_vec[i] <- .5 * se_0 ^ 2 + .5 * se_1 ^ 2
    mse_additive_vec[i] <- .5 * (
      se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    ) + .5 * (se_additive ^ 2 + (exp_additive - fx_e1) ^ 2)

    mse_e0_GxE_vec[i] <- se_0 ^ 2
    var_e0_GxE_vec[i] <- se_0 ^ 2
    bias_e0_GxE_vec[i] <- 0

    mse_e1_GxE_vec[i] <- se_1 ^ 2
    var_e1_GxE_vec[i] <- se_1 ^ 2
    bias_e1_GxE_vec[i] <- 0

    mse_e0_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    var_e0_additive_vec[i] <- se_additive ^ 2
    bias_e0_additive_vec[i] <- exp_additive - fx_e0

    mse_e1_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e1) ^ 2
    var_e1_additive_vec[i] <- se_additive ^ 2
    bias_e1_additive_vec[i] <- exp_additive - fx_e1


  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = mse_GxE_vec,
      mse_additive = mse_additive_vec,

      mse_e0_GxE = mse_e0_GxE_vec,
      var_e0_GxE = var_e0_GxE_vec,
      bias_e0_GxE = bias_e0_GxE_vec,

      mse_e1_GxE = mse_e1_GxE_vec,
      var_e1_GxE = var_e1_GxE_vec,
      bias_e1_GxE = bias_e1_GxE_vec,

      mse_e0_additive = mse_e0_additive_vec,
      var_e0_additive = var_e0_additive_vec,
      bias_e0_additive = bias_e0_additive_vec,

      mse_e1_additive = mse_e1_additive_vec,
      var_e1_additive = var_e1_additive_vec,
      bias_e1_additive = bias_e1_additive_vec,

      avg_se = avg_se_vec
    )

  return(sim_df)

}

#' Simulate MSE for GxE and Additive Model on A Grid With Difference Effect Sizes
#' And Environmental Variance Ratios
#'
#' This code allows the user to specify a grid of ratios of environmental
#' variance between environments
#' and a grid of differences in effect size between environments
#' (normalized by the standard error).
#'
#' @param sigma_ratio ratio of standard errors between environments
#' @param fx_diff_grid grid of differences in effect size
#' @param se_grid grid of standard errors of regression coefficients
#'
#' @return dataframe with various metrics for the GxE and additive models.
#' @export
#'
#' @importFrom magrittr %>%
#'
sim_var_ratio_grid_2env <- function(
  sigma_ratio_grid,
  fx_diff_grid,
  se = 1
) {

  sim_df <- expand.grid(fx_diff_grid, sigma_ratio_grid)
  colnames(sim_df) <- c("fx_diff", "sigma_ratio")

  mse_GxE_vec <- numeric(nrow(sim_df))
  mse_additive_vec <- numeric(nrow(sim_df))

  mse_e0_GxE_vec <- numeric(nrow(sim_df))
  var_e0_GxE_vec <- numeric(nrow(sim_df))
  bias_e0_GxE_vec <- numeric(nrow(sim_df))

  mse_e1_GxE_vec <- numeric(nrow(sim_df))
  var_e1_GxE_vec <- numeric(nrow(sim_df))
  bias_e1_GxE_vec <- numeric(nrow(sim_df))

  mse_e0_additive_vec <- numeric(nrow(sim_df))
  var_e0_additive_vec <- numeric(nrow(sim_df))
  bias_e0_additive_vec <- numeric(nrow(sim_df))

  mse_e1_additive_vec <- numeric(nrow(sim_df))
  var_e1_additive_vec <- numeric(nrow(sim_df))
  bias_e1_additive_vec <- numeric(nrow(sim_df))

  avg_se_vec <- numeric(nrow(sim_df))

  for (i in 1:nrow(sim_df)) {

    fx_e0 <- 0
    fx_e1 <- sim_df$fx_diff[i]

    se_0 <- se
    se_1 <- sim_df$sigma_ratio[i]

    se_additive <- sqrt(.25 * (se_0 ^ 2) + .25 * (se_1 ^ 2))

    avg_se_vec[i] <- se_additive

    exp_additive <- mean(c(fx_e0, fx_e1))

    # update vectors
    mse_GxE_vec[i] <- .5 * se_0 ^ 2 + .5 * se_1 ^ 2
    mse_additive_vec[i] <- .5 * (
      se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    ) + .5 * (se_additive ^ 2 + (exp_additive - fx_e1) ^ 2)

    mse_e0_GxE_vec[i] <- se_0 ^ 2
    var_e0_GxE_vec[i] <- se_0 ^ 2
    bias_e0_GxE_vec[i] <- 0

    mse_e1_GxE_vec[i] <- se_1 ^ 2
    var_e1_GxE_vec[i] <- se_1 ^ 2
    bias_e1_GxE_vec[i] <- 0

    mse_e0_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    var_e0_additive_vec[i] <- se_additive ^ 2
    bias_e0_additive_vec[i] <- exp_additive - fx_e0

    mse_e1_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e1) ^ 2
    var_e1_additive_vec[i] <- se_additive ^ 2
    bias_e1_additive_vec[i] <- exp_additive - fx_e1


  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = mse_GxE_vec,
      mse_additive = mse_additive_vec,

      mse_e0_GxE = mse_e0_GxE_vec,
      var_e0_GxE = var_e0_GxE_vec,
      bias_e0_GxE = bias_e0_GxE_vec,

      mse_e1_GxE = mse_e1_GxE_vec,
      var_e1_GxE = var_e1_GxE_vec,
      bias_e1_GxE = bias_e1_GxE_vec,

      mse_e0_additive = mse_e0_additive_vec,
      var_e0_additive = var_e0_additive_vec,
      bias_e0_additive = bias_e0_additive_vec,

      mse_e1_additive = mse_e1_additive_vec,
      var_e1_additive = var_e1_additive_vec,
      bias_e1_additive = bias_e1_additive_vec,

      avg_se = avg_se_vec
    )

  return(sim_df)

}

#' Simulate MSE for GxE and Additive Model on A Grid With Difference Effect Sizes
#' And Environmental Variance Ratios
#'
#' This code allows the user to specify a grid of ratios of environmental
#' variance between environments
#' and a grid of differences in effect size between environments
#' (normalized by the standard error).
#'
#' @param sigma_ratio ratio of standard errors between environments
#' @param fx_diff_grid grid of differences in effect size
#' @param se_grid grid of standard errors of regression coefficients
#'
#' @return dataframe with various metrics for the GxE and additive models.
#' @export
#'
#' @importFrom magrittr %>%
#'
get_sigma_ratio_line <- function(sigma_ratio) {

  fx_diff_grid <- seq(0, 2.5, length.out = 500)
  se_grid <- seq(1e-8, 2.5, length.out = 500)

  sim_df <- expand.grid(fx_diff_grid, se_grid)
  colnames(sim_df) <- c("fx_diff", "se")

  mse_GxE_vec <- numeric(nrow(sim_df))
  mse_additive_vec <- numeric(nrow(sim_df))

  mse_e0_GxE_vec <- numeric(nrow(sim_df))
  var_e0_GxE_vec <- numeric(nrow(sim_df))
  bias_e0_GxE_vec <- numeric(nrow(sim_df))

  mse_e1_GxE_vec <- numeric(nrow(sim_df))
  var_e1_GxE_vec <- numeric(nrow(sim_df))
  bias_e1_GxE_vec <- numeric(nrow(sim_df))

  mse_e0_additive_vec <- numeric(nrow(sim_df))
  var_e0_additive_vec <- numeric(nrow(sim_df))
  bias_e0_additive_vec <- numeric(nrow(sim_df))

  mse_e1_additive_vec <- numeric(nrow(sim_df))
  var_e1_additive_vec <- numeric(nrow(sim_df))
  bias_e1_additive_vec <- numeric(nrow(sim_df))

  avg_se_vec <- numeric(nrow(sim_df))

  for (i in 1:nrow(sim_df)) {

    fx_e0 <- 0
    fx_e1 <- sim_df$fx_diff[i]

    se_0 <- sim_df$se[i]
    se_1 <- se_0 * sigma_ratio

    se_additive <- sqrt(.25 * (se_0 ^ 2) + .25 * (se_1 ^ 2))

    avg_se_vec[i] <- se_additive

    exp_additive <- mean(c(fx_e0, fx_e1))

    # update vectors
    mse_GxE_vec[i] <- .5 * se_0 ^ 2 + .5 * se_1 ^ 2
    mse_additive_vec[i] <- .5 * (
      se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    ) + .5 * (se_additive ^ 2 + (exp_additive - fx_e1) ^ 2)

    mse_e0_GxE_vec[i] <- se_0 ^ 2
    var_e0_GxE_vec[i] <- se_0 ^ 2
    bias_e0_GxE_vec[i] <- 0

    mse_e1_GxE_vec[i] <- se_1 ^ 2
    var_e1_GxE_vec[i] <- se_1 ^ 2
    bias_e1_GxE_vec[i] <- 0

    mse_e0_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e0) ^ 2
    var_e0_additive_vec[i] <- se_additive ^ 2
    bias_e0_additive_vec[i] <- exp_additive - fx_e0

    mse_e1_additive_vec[i] <- se_additive ^ 2 + (exp_additive - fx_e1) ^ 2
    var_e1_additive_vec[i] <- se_additive ^ 2
    bias_e1_additive_vec[i] <- exp_additive - fx_e1


  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = mse_GxE_vec,
      mse_additive = mse_additive_vec,

      mse_e0_GxE = mse_e0_GxE_vec,
      var_e0_GxE = var_e0_GxE_vec,
      bias_e0_GxE = bias_e0_GxE_vec,

      mse_e1_GxE = mse_e1_GxE_vec,
      var_e1_GxE = var_e1_GxE_vec,
      bias_e1_GxE = bias_e1_GxE_vec,

      mse_e0_additive = mse_e0_additive_vec,
      var_e0_additive = var_e0_additive_vec,
      bias_e0_additive = bias_e0_additive_vec,

      mse_e1_additive = mse_e1_additive_vec,
      var_e1_additive = var_e1_additive_vec,
      bias_e1_additive = bias_e1_additive_vec,

      avg_se = avg_se_vec
    )

  sim_df <- sim_df %>%
    dplyr::mutate(abs_mse_diff = abs(mse_e0_additive - mse_e0_GxE))

  quant <- quantile(sim_df$abs_mse_diff, .005)

  reg_df <- sim_df %>%
    dplyr::filter(abs_mse_diff <= quant)

  mod <- lm(se ~ fx_diff, data = reg_df)

  coefs <- coef(summary(mod))

  return(
    list(
      intercept = coefs["(Intercept)", "Estimate"],
      slope = coefs["fx_diff", "Estimate"]
    )
  )


}

# Here, want to build a simulation where I can set a true difference in effects
# and then I can vary the noise
# as I vary the noise, we will see that as the estimated difference in effects is larger
# the decision rule will shift
sim_beta_hat_diff <- function(diff_range, sd_range, pts_per_sim) {
  
  '%>%' <- magrittr::'%>%'
  sim_df <- expand.grid(diff_range, sd_range)
  colnames(sim_df) <- c("diff", "sd")
  sim_df <- sim_df %>%
    dplyr::mutate(e0_fx = 0, e1_fx = diff)
  
  sim_df <- do.call("rbind", replicate(pts_per_sim, sim_df, simplify = FALSE))
  
  e0_mles <- rnorm(
    n = nrow(sim_df), 
    sd = sim_df$sd,
    mean = sim_df$e0_fx
  )
  
  e1_mles <- rnorm(
    n = nrow(sim_df),
    sd = sim_df$sd,
    mean = sim_df$e1_fx
  )
  
  add_mles <- .5 * (e0_mles + e1_mles)
  
  sim_df <- sim_df %>% 
    dplyr::mutate(
      e0_est = e0_mles,
      e1_est = e1_mles,
      add_est = add_mles
  )
  
  ash_subset <- sim_df %>%
    dplyr::select(
      e0_est, e1_est, add_est, sd
    ) %>%
    dplyr::sample_n(min(100000, nrow(sim_df)))

  ash_model <- ashr::ash(
    betahat = ash_subset$e1_est - ash_subset$e0_est,
    sebetahat = sqrt(2) * ash_subset$sd,
    mixcompdist = "normal"
  )

  g_fit <- ashr::get_fitted_g(ash_model)
  g_fit$pi <- rep(0, length(g_fit$pi))
  g_fit$pi[1] <- .5
  g_fit$pi[4] <- .25
  g_fit$pi[17] <- .25
  
  ash_posteriors <- ashr::ash(
    betahat = sim_df$e1_est - sim_df$e0_est,
    sebetahat = sqrt(2) * sim_df$sd,
    fixg = TRUE,
    g = g_fit,
    outputlevel = 1
  )
  
  sim_df <- sim_df %>%
    dplyr::mutate(
      e0_GxE_error = (e0_est - e0_fx) ^ 2,
      e1_GxE_error = (e1_est - e1_fx) ^ 2,
      e0_add_error = (add_est - e0_fx) ^ 2,
      e1_add_error = (add_est - e1_fx) ^ 2,
      abs_beta_hat_diff = abs(e1_est - e0_est),
      abs_beta_hat_diff_bc = pmax(0, abs(e1_est - e0_est) - sqrt(4 / pi) * sd),
      ash_beta_diff = abs(ashr::get_pm(ash_posteriors)),
      ash_beta_diff_sd = ashr::get_psd(ash_posteriors)
    )
  
  ash_diff_sim_df <- sim_df %>%
    dplyr::filter(
      ash_beta_diff <= .3
    )
  
  # exclude all points with a beta hat diff more than 2
  beta_hat_diff_sim_df <- sim_df %>%
    dplyr::filter(abs_beta_hat_diff <= 1.5)
  
  beta_hat_bc_sim_df <- sim_df %>%
    dplyr::filter(abs_beta_hat_diff_bc < 1)
  
  beta_hat_diff_sim_df <- beta_hat_diff_sim_df %>%
    dplyr::mutate(
      percentile = cut(
        abs_beta_hat_diff,
        seq(0, 1.6, length.out = 100),
        include.lowest = TRUE,
        labels = FALSE
      )
    )
  
  beta_hat_bc_sim_df <- beta_hat_bc_sim_df %>%
    dplyr::mutate(
      percentile = cut(
        abs_beta_hat_diff_bc,
        seq(0, 1, length.out = 100),
        include.lowest = TRUE,
        labels = FALSE
      )
    )
  
  # ash_diff_sim_df <- ash_diff_sim_df %>%
  #   dplyr::mutate(
  #     percentile = cut(
  #       ash_beta_diff,
  #       seq(0, .3, length.out = 50),
  #       include.lowest = TRUE,
  #       labels = FALSE
  #     )
  #   )
  
  # ash_diff_sim_df <- ash_diff_sim_df %>%
  #   dplyr::mutate(
  #     percentile_sd = cut(
  #       ash_beta_diff_sd,
  #       seq(0.25, 0.55, length.out = 50),
  #       include.lowest = TRUE,
  #       labels = FALSE
  #     )
  #   )
  
  beta_hat_diff_sim_df <- beta_hat_diff_sim_df %>%
    dplyr::mutate(
      mean_abs_diff = seq(0, 1.5, length.out = 100)[percentile]
    )
  
  beta_hat_bc_sim_df <- beta_hat_bc_sim_df %>%
    dplyr::mutate(
      mean_abs_bc_diff = seq(0, 1, length.out = 100)[percentile]
    )
  
  # ash_diff_sim_df <- ash_diff_sim_df %>%
  #   dplyr::mutate(
  #     mean_ash_diff = seq(0, .3, length.out = 50)[percentile]
  #     #mean_ash_sd_diff = seq(0.25, 0.55, length.out = 50)[percentile_sd]
  #   )
  
  beta_hat_diff_sum_sim_df <- beta_hat_diff_sim_df %>%
    dplyr::group_by(sd, mean_abs_diff) %>%
    dplyr::summarize(
      e0_GxE_mse = mean(e0_GxE_error),
      e1_GxE_mse = mean(e1_GxE_error),
      e0_add_mse = mean(e0_add_error),
      e1_add_mse = mean(e1_add_error),
      sample_size = dplyr::n()
    )
  
  beta_hat_bc_sum_sim_df <- beta_hat_bc_sim_df %>%
    dplyr::group_by(sd, mean_abs_bc_diff) %>%
    dplyr::summarize(
      e0_GxE_mse = mean(e0_GxE_error),
      e1_GxE_mse = mean(e1_GxE_error),
      e0_add_mse = mean(e0_add_error),
      e1_add_mse = mean(e1_add_error),
      sample_size = dplyr::n()
    )
  
  # ash_diff_sum_sim_df <- ash_diff_sim_df %>%
  #   dplyr::group_by(sd, mean_ash_diff) %>%
  #   dplyr::summarize(
  #     e0_GxE_mse = mean(e0_GxE_error),
  #     e1_GxE_mse = mean(e1_GxE_error),
  #     e0_add_mse = mean(e0_add_error),
  #     e1_add_mse = mean(e1_add_error),
  #     sample_size = dplyr::n()
  #   )
  
  return(
    list(
      beta_hat_diff_df = beta_hat_diff_sum_sim_df,
      beta_hat_bc_sum_sim_df = beta_hat_bc_sum_sim_df
    )
   )
  
}

