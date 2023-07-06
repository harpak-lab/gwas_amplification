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
    se_1 <- se_0 / sigma_ratio

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

sim_var_ratio_2env_v2 <- function(
    r,
    fx_diff_grid,
    se_grid
) {
  
  sim_df <- expand.grid(fx_diff_grid, se_grid)
  colnames(sim_df) <- c("fx_diff", "se_A")
  
  sim_df <- sim_df %>%
    dplyr::mutate(V_A = se_A ^ 2)
  
  sim_df <- sim_df %>%
    dplyr::mutate(V_B = V_A / r)
  
  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_GxE = V_A,
      mse_additive = .25 * V_A + .25 * V_B + (fx_diff / 2) ^ 2
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
get_sigma_ratio_slope <- function(r, n, m) {

  slope <- m / (
    (n + m) * sqrt(1 - ((n / (n + m)) ^ 2) - ((m ^ 2) / (r * ((n + m) ^ 2))))
  )
  
  return(slope)

}

make_cv_into_grid <- function(
    cv_df, 
    se_start, 
    se_end, 
    se_step,
    fx_diff_start,
    fx_diff_end,
    fx_diff_step
  ) {
  
  fx_diff_grid <- seq(fx_diff_start + fx_diff_step, fx_diff_end, fx_diff_step)
  se_grid <- seq(se_start + se_step, se_end, se_step)
  sim_df <- expand.grid(fx_diff_grid, se_grid)
  colnames(sim_df) <- c("fx_diff", "se")
  
  mse_diff_vec <- numeric(nrow(sim_df))
  samp_size_vec <- numeric(nrow(sim_df))
  
  # now, I know exactly how much to step in each direction
  ctr <- 1
  
  for (i in 1:nrow(sim_df)) {
    
    print(glue::glue("performing {i} of {nrow(sim_df)} total operations"))
    grid_pt_df <- cv_df %>%
      dplyr::filter(
        se < sim_df$se[i] & se >= (sim_df$se[i] - se_step) &
          fx_diff < sim_df$fx_diff[i] & fx_diff >= (sim_df$fx_diff[i] - fx_diff_step)
      )
    
    if (nrow(grid_pt_df) == 0) {
      
      mse_diff_vec[i] <- 0
      samp_size_vec[i] <- 0
      
    } else {
      
      samp_size_vec[i] <- nrow(grid_pt_df)
      mse_diff_vec[i] <- mean(grid_pt_df$mse_diff)
      
    }
    
  }
  
  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff = mse_diff_vec,
      sample_size = samp_size_vec
    )
  
  return(sim_df)
  
}
