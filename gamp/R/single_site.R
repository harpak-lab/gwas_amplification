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
