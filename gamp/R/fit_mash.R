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
generate_mash_fx <- function(n, p, Sigma, pi, s, mu = NULL) {

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

  X <- matrix(0, n, p+1)
  for (i in nonnull_fx) {

    if(is.null(mu)) {
      
      X[i,] <- c(
        MASS::mvrnorm(
          1, rep(0, p), Sigma[[which_sigma[i]]]), which_sigma[i]
        )
      
    } else {
      
      X[i,] <- c(
        MASS::mvrnorm(
          1, mu[[which_sigma[i]]], Sigma[[which_sigma[i]]]), which_sigma[i]
      )
      
    }

  }

  return(X)

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


#' Fit mash to 2 Environment GWAS Data
#'
#' mash is fit with canonical covariance matrices only, where each matrix is one
#' element of the cartesian produce of \code{corr_range} and \code{amp_range}.
#' See \code{make_amp_cov_mat_list} for more details.
#'
#' @param Bhat n x 2 matrix of estimated effect sizes
#' @param Shat n x 2 matrix of standard errors
#' @param corr_range vector indicating range of correlation for covariance matrices.
#' Must be beteween 0 and 1.
#' @param amp_range vector indicating range of amplification of covariance matrices.
#'
#' @return mash object
#'
fit_mash_gwas_2env <- function(
    Bhat, Shat, cov_mats = NULL, corr_range = NULL, amp_range = NULL, fitted_g = NULL
  ) {

  if (is.null(cov_mats)) {
    
    cov_mats <- make_amp_cov_mat_list(corr_range, amp_range)
    
  }

  mash_data <- mashr::mash_set_data(Bhat, Shat)
  if (is.null(fitted_g)) {

    m <- mashr::mash(mash_data, cov_mats)

  } else {

    m <- mashr::mash(mash_data, fixg = TRUE, g = fitted_g)

  }

  return(m)

}
