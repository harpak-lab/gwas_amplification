get_ssg_df <- function(n, maf_quantiles, num_sims = 5) {

  ssg_vec <- c()

  for (i in 2:length(maf_quantiles)) {

    ssg_est <- 0
    reg_coef_est <- 0
    maf <- maf_quantiles[i - 1] * .5 + maf_quantiles[i] * .5

    for (j in 1:num_sims) {

      g <- rbinom(n, 2, maf)

      g_bar <- mean(g)
      ssg_est <- ssg_est + (1 / num_sims) * sum((g - g_bar) ^ 2)

    }

    ssg_vec <- c(ssg_vec, ssg_est)

  }

  sim_df <- data.frame(
    maf_low = maf_quantiles[1:(length(maf_quantiles) - 1)],
    maf_high = maf_quantiles[2:length(maf_quantiles)],
    ssg = ssg_vec
  )

  return(sim_df)

}

