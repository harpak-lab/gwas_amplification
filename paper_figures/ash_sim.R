sim_ash_train_test_dif <- function(
  train_dist, test_dist, n_train, n_test, n_sims
) {

  test_mse_vec <- c()

  for (i in 1:n_sims) {

    train_sampler <- distr::r(train_dist)
    train_samples <- train_sampler(n_train) + rnorm(n_train)

    test_sampler <- distr::r(test_dist)
    test_x <- test_sampler(n_test)
    test_samples <- test_x + rnorm(n_test)

    ash_train <- ashr::ash(
      train_samples, rep(1, n_train), mixcompdist = "normal"
    )

    ash_test <- ashr::ash(
      test_samples,
      rep(1, n_train),
      fixg = TRUE,
      g = ashr::get_fitted_g(ash_train),
      mixcompdist = "normal"
    )

    ash_preds <- ash_test$result$PosteriorMean
    test_mse <- mean((ash_preds - test_x) ^ 2)
    test_mse_vec <- c(test_mse_vec, test_mse)

  }

  return(mean(test_mse_vec))

}

mix1 <- distr::UnivarMixingDistribution(
  distr::Norm(mean = 0, sd = 1),
  distr::Dirac(),
  mixCoeff = c(.35, .65)
)

sim1_out <- sim_ash_train_test_dif(
  mix1, mix1, n_train = 10000, n_test = 10000, n_sims = 10
)

mix2 <- distr::UnivarMixingDistribution(
  distr::Norm(mean = 0, sd = 1),
  distr::Dirac(),
  mixCoeff = c(.1, .9)
)

sim2_out <- sim_ash_train_test_dif(
  mix2, mix1, n_train = 10000, n_test = 10000, n_sims = 10
)
