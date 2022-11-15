# code for figure 3 - analysis of pgs as % amplification increases
set.seed(1)
'%>%' <- magrittr::'%>%'
#library(ggplot2)
#library(foreach)

load("~/Downloads/mash_100_R/diastolicBP_auto_mash_100g.RData")

get_Sigma_from_fitted_g <- function(g) {

  Sigma <- list()
  num_components <- length(g$pi)

  mixture_comp_ctr <- 1

  # sparsity
  s <- g$pi[1]

  pi_out <- c()

  for (i in 2:num_components) {

    if (g$pi[i] > 0) {

      tokenized <- stringr::str_split(names(g_ave$pi)[i], pattern = "\\.")[[1]]
      grid_pt <- as.numeric(tokenized[length(tokenized)])
      cov_mat_name <- paste(tokenized[-length(tokenized)], collapse = ".")
      cov_mat <- g$Ulist[[cov_mat_name]]
      Sigma[[mixture_comp_ctr]] <- g$grid[grid_pt] * cov_mat
      pi_out <- c(pi_out, g$pi[i])
      mixture_comp_ctr <- mixture_comp_ctr + 1

    }

  }

  return(
    list(
      Sigma = Sigma,
      s = s,
      pi = pi_out
    )
  )

}

out <- get_Sigma_from_fitted_g(g_list[[5]])

generate_fig5 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {


    maf_simulator <- function(n) {

      pmin(rbeta(n, 1, 10), .5 - 1e-05)

    }

    sim_out <- gamp::simulate_pgs_corr(
      n_e0 = 1.5e05,
      n_e1 = 1.5e05,
      n_snps = 1700,
      maf_simulator = maf_simulator,
      e0_h2 = .25,
      e1_h2 = .25,
      Sigma = out$Sigma,
      pi = out$pi,
      s = out$s,
      num_sims = 10
    )

  }

  if(write_rds) {

    readr::write_rds(sim_out, "./rds_data/fig5_bp.rds")

  }

  # now, convert the dataframe to long format for ggplot
  # implement once I actually have the dataframe on hand


}

# first time simulation
generate_fig5(write_rds = TRUE)
