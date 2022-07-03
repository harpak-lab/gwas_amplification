# sample genotypes given number of individuals and minor allele frequency
sample_genotypes <- function(n, maf) {

  sample(
    x=c(0, 1, 2),
    size=n,
    prob=c((1 - maf) ^ 2, 2 * maf * (1 - maf), maf ^ 2),
    replace = TRUE
  )

}

#' Generate Data from Mash Model
#'
#' @param n number of samples
#' @param p number of conditions
#' @param Sigma list of covariance matrices
#' @param pi probability of sampling from each covariance matrix
#' @param s percent of the signals to be 0 for all conditions
#'
#' @return Matrix of true effects X
#' @export
#'
#' @examples
generate_mash_fx <- function(n, p, Sigma, pi, s) {

  if (is.null(pi)) {
    pi = rep(1, length(Sigma)) # default to uniform distribution
  }
  assertthat::are_equal(length(pi), length(Sigma))
  for (j in length(Sigma)) {
    assertthat::are_equal(dim(Sigma[j]), c(p, p))
  }

  pi <- pi / sum(pi) # normalize pi to sum to one
  which_sigma <- sample(1:length(pi), n, replace=TRUE, prob=pi)
  nonnull_fx <- sample(1:n, floor((1 - s)*n), replace=FALSE)

  X <- matrix(0, n, p)
  for (i in nonnull_fx) {
    X[i,] <- MASS::mvrnorm(1, rep(0, p), Sigma[[which_sigma[i]]])
  }

  return(X)

}

# helper function to generate dataframe in two environments
simulate_two_env_df <- function(
  n_e0,
  n_e1,
  fx_e0,
  fx_e1,
  sigma_e0,
  sigma_e1,
  int_e0,
  int_e1,
  maf_e0,
  maf_e1
) {

  g_e0 <- sample_genotypes(n_e0, maf_e0)

  g_e1 <- sample_genotypes(n_e1, maf_e1)

  y_e0 <- int_e0 + g_e0 * fx_e0 + rnorm(n = n_e0, sd = sigma_e0)
  y_e1 <- int_e1 + g_e1 * fx_e1 + rnorm(n = n_e1, sd = sigma_e1)

  sim_df <- data.frame(
    y = c(y_e0, y_e1),
    g = c(g_e0, g_e1),
    e = c(rep("e0", n_e0), rep("e1", n_e1))
  )

  return(sim_df)

}

#' Generate GWAS Summary Statistics For Two Correlated Environments
#'
#' @param X m x 2 matrix of effect sizes where there are m snps
#' @param n_e0 number of individuals in environment 0
#' @param n_e1 number of individuals in environment 1
#' @param int_e0 intercept for environment 0
#' @param int_e1 intercept for environment 1
#' @param maf_e0 minor allele frequency in environment 0
#' @param maf_e1 minor allele frequency in environment 1
#' @param sigma_e0 observation noise in environment 0
#' @param sigma_e1 observation noise in environment 1
#'
#' @return list of two matrices. Bhat are effect sizes and Shat are standard
#' errors
#' @export
#'
#' @examples
generate_2env_gwas <- function(
  X,
  n_e0,
  n_e1,
  sigma_e0,
  sigma_e1,
  int_e0,
  int_e1,
  maf_e0,
  maf_e1
) {

  m <- nrow(X)

  Bhat_sep <- matrix(nrow = m, ncol = 2)
  Shat_sep <- matrix(nrow = m, ncol = 2)

  Bhat_comb <- matrix(nrow = m, ncol = 1)
  Shat_comb <- matrix(nrow = m, ncol = 1)


  for (i in 1:m) {

    snp_df <- simulate_two_env_df(
      n_e0 = n_e0,
      n_e1 = n_e1,
      sigma_e0 = sigma_e0,
      sigma_e1 = sigma_e1,
      int_e0 = int_e0,
      int_e1 = int_e1,
      maf_e0 = maf_e0,
      maf_e1 = maf_e1,
      fx_e0 = X[i, 1],
      fx_e1 = X[i, 2]
    )

    comb_lm <- lm(y ~ g, data = snp_df)
    Bhat_comb[i, 1] <- coef(summary(comb_lm))['g', 'Estimate']
    Shat_comb[i, 1] <- coef(summary(comb_lm))['g', 'Std. Error']

    e0_lm <- lm(y ~ g, data = snp_df, subset = c(e == "e0"))
    e1_lm <- lm(y ~ g, data = snp_df, subset = c(e == "e1"))

    Bhat_sep[i, 1] <- coef(summary(e0_lm))['g', 'Estimate']
    Shat_sep[i, 1] <- coef(summary(e0_lm))['g', 'Std. Error']

    Bhat_sep[i, 2] <- coef(summary(e1_lm))['g', 'Estimate']
    Shat_sep[i, 2] <- coef(summary(e1_lm))['g', 'Std. Error']

  }

  return(
    list(
      Bhat_comb = Bhat_comb,
      Shat_comb = Shat_comb,
      Bhat_sep = Bhat_sep,
      Shat_sep = Shat_sep
    )
  )

}

#' Create 2x2 covariance matrices with optional levels of amplification
#'
#' The function assumes that we have two groups, with the standard deviation of
#' the effect size in one group being 1, and the standard deviation of the effect
#' size in the other group being \code{1 * amp_coef} if \code{amp} is set to
#' \code{TRUE}
#'
#' @param desired_corr Desired level of correlation. Must be in [-1, 1].
#' @param amp_coef Coefficient of amplification, as described above.
#' @param amp Boolean indicating if any amplificaiton should take place
#' @param amp_hs Boolean indicating if amplification should take place in group
#' e0 or group e1. Only used if \code{amp} is set to \code{TRUE}.
#'
#' @return 2x2 covariance matrix
#' @export
#'
#' @examples
make_amp_cov_mat <- function(
  desired_corr, amp_coef = 1, amp = TRUE, amp_e0 = TRUE
) {

  if (amp_e0 && amp) {

    e1_sd <- 1
    e0_sd <- e1_sd * amp_coef

  } else if(!amp_e0 && amp) {

    e0_sd <- 1
    e1_sd <- e0_sd * amp_coef

  } else {

    e0_sd <- 1
    e1_sd <- 1

  }

  # derive covariance from correlation and sds
  cov_e0_e1 <- desired_corr * e0_sd * e1_sd

  cov_mat <- matrix(
    data = c(e0_sd ^ 2, cov_e0_e1, cov_e0_e1, e1_sd ^ 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      rows = c("e0", "e1"), cols = c("e0", "e1")
    )
  )

  return(cov_mat)

}

#' Create grid of amplification covariance matrices
#'
#' @param corr_range vector of discrete values of correlation between
#' environmental effects
#' @param amp_range vector of discreate values of amplification of one
#' environment over the other
#'
#' @return list of covariance matrices
#' @export
#'
#' @examples
make_amp_cov_mat_list <- function(
  corr_range,
  amp_range
) {

  cov_mat_list <- list()

  cov_mat_list[['e1_spec']] <- matrix(
    data = c(0, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(
      rows = c("e0", "e1"), cols = c("e0", "e1")
    )
  )

  cov_mat_list[['e0_spec']] <- matrix(
    data = c(1, 0, 0, 0), nrow = 2, byrow = TRUE, dimnames = list(
      rows = c("e0", "e1"), cols = c("e0", "e1")
    )
  )

  for(corr in corr_range) {

    cov_mat_list[[glue::glue('equal_corr_{corr}')]] <- make_amp_cov_mat(
      desired_corr = corr, amp = FALSE
    )

    for(cond in c("e0", "e1")) {

      for(amp in amp_range) {

        cov_mat_list[[glue::glue('{cond}_amp_{amp}_corr_{corr}')]] <- make_amp_cov_mat(
          desired_corr = corr, amp_e0 = (cond == "e0"), amp_coef = amp
        )

      }

    }

  }

  return(cov_mat_list)

}

# fit mash to gwas data
fit_mash_gwas <- function(Bhat, Shat, corr_range, amp_range) {

  cov_mats <- make_amp_cov_mat_list(corr_range, amp_range)
  mash_data <- mashr::mash_set_data(Bhat, Shat)
  m <- mashr::mash(mash_data, cov_mats)
  return(m)

}

mash_mse <- function(m, X) {
  mean((ashr::get_pm(m) - X) ^ 2)
}

gwas_mse_2env <- function(Bhat, X) {

  if (ncol(Bhat) == 2) {

    mse <- mean((Bhat - X) ^ 2)

  } else {

    Bhat_rep <- matrix(data = rep(as.vector(Bhat), 2), ncol = 2)
    mse <- mean((Bhat_rep - X) ^ 2)

  }

  return(mse)

}

generate_2env_fx_ecdf <- function(n, fx_vec, e1_amp = 1) {

  fx_samp_e0 <- sample(x = fx_vec, size = n, replace = TRUE)
  fx_samp_e1 <- e1_amp * fx_samp_e0
  fx_mat <- matrix(data = c(fx_samp_e0, fx_samp_e1), nrow = n, ncol = 2)
  return(fx_mat)

}

#' Simulate MSE of Different Modelling Techniques for Two Environment GWAS
#'
#' @param n_sims number of simulations to run
#' @param n_snps number of snps in each simulation
#' @param n_e0 number of individuals in environment 0
#' @param n_e1 number of individuals in environment 1
#' @param sigma_e0 observation standard error in environment 0
#' @param sigma_e1 observation standard error in environment 1
#' @param Sigma list of covariance matrices to generate true effects
#' @param pi (optional) vector of probabilities to determine weights on
#' matrices in Sigma. Defaults to uniform.
#' @param s (optional) level of sparsity in true effects
#' @param corr_range (optional) vector of correlation between effects in
#' environment 1 and environment 0 to consider in mash model
#' @param amp_range (optional) vector of amplification level between effects in
#' environment 1 and environment 0 to consider in mash model
#' @param int_e0 (optional) intercept in environment 0. Applies to all snps.
#' @param int_e1 (optional) intercept in environment 1. Applies to all snps.
#' @param maf_e0 (optional) minor allele frequency in environment 0.
#' Applies to all snps.
#' @param maf_e1 (optional) minor allele frequency in environment 1.
#' Applies to all snps.
#'
#' @return list of mse values for different modelling techniques.
#' @export
#'
#' @examples
simulate_2env_gwas_mse <- function(
  n_sims,
  n_snps,
  n_e0,
  n_e1,
  sigma_e0,
  sigma_e1,
  Sigma = NULL,
  prior = c("mash", "ecdf"),
  fx_ecdf = NULL,
  ecdf_e1_amp = 1,
  pi = NULL,
  s = 0.8,
  corr_range = seq(from = -1, to = 1, by = .25),
  amp_range = c(1.5, 2, 3),
  int_e0 = 0,
  int_e1 = 0,
  maf_e0 = .4,
  maf_e1 = .4
) {

  prior <- match.arg(prior)
  sep_mse_vec <- numeric(n_sims)
  comb_mse_vec <- numeric(n_sims)
  mash_mse_vec <- numeric(n_sims)

  for (i in 1:n_sims) {

    print(glue::glue("Running simulation {i} of {n_sims}..."))

    if (prior == "mash") {

      # generate true effects from mash model
      fx_mat <- generate_mash_fx(
        n = n_snps, p = 2, Sigma = Sigma, pi = pi, s = s
      )

    } else if (prior == "ecdf") {

      fx_mat <- generate_2env_fx_ecdf(n_snps, fx_ecdf, ecdf_e1_amp)

    }

    # generate gwas summary statistics from true effects
    gwas_stats <- generate_2env_gwas(
      fx_mat,
      n_e0,
      n_e1,
      sigma_e0,
      sigma_e1,
      int_e0,
      int_e1,
      maf_e0,
      maf_e1
    )

    # fit mash model to two environment gwas data
    mash_fit <- fit_mash_gwas(
      Bhat = gwas_stats$Bhat_sep,
      Shat = gwas_stats$Shat_sep,
      corr_range = corr_range,
      amp_range = amp_range
    )

    # get errors
    mash_mse_vec[i] <- mash_mse(mash_fit, fx_mat)
    comb_mse_vec[i] <- gwas_mse_2env(gwas_stats$Bhat_comb, fx_mat)
    sep_mse_vec[i] <- gwas_mse_2env(gwas_stats$Bhat_sep, fx_mat)

  }

  return(
    list(
      mash_mse = mean(mash_mse_vec),
      sep_mse = mean(sep_mse_vec),
      comb_mse = mean(comb_mse_vec)
    )
  )

}
