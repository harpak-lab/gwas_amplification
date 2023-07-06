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

protected_cor <- function(x, y) {

  if (all(x == 0)) {

    return(0)

  } else {

    return(cor(x, y))

  }

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
  n_e0_train,
  n_e1_train,
  n_e0_test,
  n_e1_test,
  n_snps,
  maf_simulator,
  e0_h2,
  e1_h2,
  Sigma,
  pi,
  n_ascertained_snps,
  num_sims = 10,
  s = 0,
  beta_methods = c("additive", "GxE", "mash", "ash_additive", "ash_GxE"),
  selection_methods = c("additive", "GxE", "all")
) {

  sim_results <- list()

  for (selection_method in selection_methods) {

    sim_results[[selection_method]] <- list()

    for (beta_method in beta_methods) {

      sim_results[[selection_method]][[beta_method]] <- list()
      sim_results[[selection_method]][[beta_method]][["corr"]] <- 0

    }

  }

  for (i in 1:num_sims) {

    # sample from mash prior
    fx_mat <- generate_mash_fx(
      n = n_snps, p = 2, pi = pi, Sigma = Sigma, s = s
    )

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

    GxE_models_e0 <- estimate_GxE_model_fast(
      n = n_e0_train, 
      maf_vec = maf_vec, 
      fx_vec = fx_mat[, 1], 
      h2 = e0_h2
    )
    
    GxE_models_e1 <- estimate_GxE_model_fast(
      n = n_e1_train, 
      maf_vec = maf_vec, 
      fx_vec = fx_mat[, 2], 
      h2 = e1_h2
    )
    
    GxE_models <- list(
      est_e0 = GxE_models_e0$bhat,
      sd_e0 = GxE_models_e0$shat,
      pval_e0 = GxE_models_e0$pval,
      est_e1 = GxE_models_e1$bhat,
      sd_e1 = GxE_models_e1$shat,
      pval_e1 = GxE_models_e1$pval
    )
    
    additive_models <- list(
      est = .5 * GxE_models$est_e0 + .5 * GxE_models$est_e1,
      sd = sqrt(.25 * (GxE_models$sd_e0 ^ 2) + .25 * (GxE_models$sd_e1 ^ 2))
    )
    
    additive_models[['pval']] <- pt(
      -abs(additive_models$est / additive_models$sd), df = n_e0_train + n_e1_train - 2
      ) + 
      (1 - pt(abs(additive_models$est / additive_models$sd), df = n_e0_train + n_e1_train - 2))

    selected_snps <- list()
    beta_ests <- list()
    
    sorted_pvals_GxE_e0 <- sort(GxE_models$pval_e0)
    sorted_pvals_GxE_e1 <- sort(GxE_models$pval_e1)
    
    pval_thresh_GxE_e0 <- sorted_pvals_GxE_e0[n_ascertained_snps]
    pval_thresh_GxE_e1 <- sorted_pvals_GxE_e1[n_ascertained_snps]
    
    selected_snps[["GxE"]] <- matrix(
      data = c(
        GxE_models$pval_e0 <= pval_thresh_GxE_e0,
        GxE_models$pval_e1 <= pval_thresh_GxE_e1
      ),
      ncol = 2
    )
    
    sorted_pvals_additive <- sort(additive_models$pval)
    pval_thresh_additive <- sorted_pvals_additive[n_ascertained_snps]
    
    selected_snps[["additive"]] <- matrix(
      data = rep(additive_models$pval <= pval_thresh_additive, 2), ncol = 2
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

    if ("mash" %in% beta_methods || "mash" %in% selection_methods) {

      skip_to_next <- FALSE
      
      tryCatch(
        {
          
          mash_model <- get_test_preds_mash(
            Bhat, Shat, test_df, amp_range = c(1.5, 2, 3)
          )
          
        }, error = function(cond) {
          
          skip_to_next <<- TRUE
          
        }
      )
      
      if(skip_to_next) {
        
        next
        
      }
      
      beta_ests[['mash']] <- matrix(
        data = c(mash_model$est_e0, mash_model$est_e1),
        ncol = 2
      )
      
      sorted_lfsr_e0 <- sort(mash_model$lfsr_e0)
      sorted_lfsr_e1 <- sort(mash_model$lfsr_e1)
      
      lfsr_thresh_e0 <- sorted_lfsr_e0[n_ascertained_snps]
      lfsr_thresh_e1 <- sorted_lfsr_e1[n_ascertained_snps]
      
      selected_snps[["mash"]] <- matrix(
        data = c(
          mash_model$lfsr_e0 <= lfsr_thresh_e0,
          mash_model$lfsr_e1 <= lfsr_thresh_e1
        ),
        ncol = 2
      )

    }
    
    if ("ash_additive" %in% beta_methods) {
      
      skip_to_next <- FALSE
      
      tryCatch(
        
        {
          
          ash_model_add <- ashr::ash(additive_models$est, additive_models$sd)
          post_ests_add <- ashr::get_pm(ash_model_add)
          
        }, error = function(cond) {
          
          skip_to_next <<- TRUE
          
        }
        
      )
      
      if(skip_to_next) {
        
        next
        
      }
      
      beta_ests[['ash_additive']] <- matrix(
        data = rep(post_ests_add, 2), ncol = 2
      )
      
    }
    
    if ("ash_GxE" %in% beta_methods) {
      
      skip_to_next <- FALSE
      
      tryCatch(
        {
          
          ash_model_GxE_e0 <- ashr::ash(GxE_models$est_e0, GxE_models$sd_e0)
          post_ests_GxE_e0 <- ashr::get_pm(ash_model_GxE_e0)
          
          ash_model_GxE_e1 <- ashr::ash(GxE_models$est_e1, GxE_models$sd_e1)
          post_ests_GxE_e1 <- ashr::get_pm(ash_model_GxE_e1)
          
        }, error = function(cond) {
          
          skip_to_next <<- TRUE
          
        }
      )
      
      if(skip_to_next) {
        
        next
        
      }
      
      beta_ests[['ash_GxE']] <- matrix(
        data = c(post_ests_GxE_e0, post_ests_GxE_e1), ncol = 2
      )
      
    }

    for (selection_method in selection_methods) {
      
      for (beta_method in beta_methods) {
        
        # get the correlation for each of the prediction methods
        selected_beta_ests <- beta_ests[[beta_method]]
        
        if (selection_method != "all") {
          
          selected_beta_ests[!selected_snps[[selection_method]]] <- 0
          
        }

        preds <- c(
          test_mat_e0 %*% selected_beta_ests[, 1],
          test_mat_e1 %*% selected_beta_ests[, 2]
        )
        
        sim_results[[selection_method]][[beta_method]][['corr']] <- sim_results[[selection_method]][[beta_method]][['corr']] +
          (1 / num_sims) * protected_cor(preds, test_df$y)
        
      }
      
    }

  }

  return(sim_results)

}
