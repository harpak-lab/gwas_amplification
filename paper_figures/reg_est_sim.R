# I want to roughly verify my theoretical calculation via simulation

sim_add_reg_est <- function(
  n,
  m,
  f_A,
  f_B, 
  beta_A,
  beta_B,
  sigma_A,
  sigma_B
  ) {
  
  g_A <- rbinom(n = n, size = 2, prob = f_A)
  g_A_bar <- mean(g_A)
  g_B <- rbinom(n = m, size = 2, prob = f_B)
  g_B_bar <- mean(g_B)
  
  y_A <- g_A * beta_A + rnorm(n = n, sd = sigma_A)
  y_B <- g_B * beta_B + rnorm(n = m, sd = sigma_B)
  
  y_A <- scale(y_A, scale = FALSE)
  y_B <- scale(y_B, scale = FALSE)
  
  lm_A <- lm(y_A ~ g_A)
  beta_hat_A <- coef(summary(lm_A))["g_A", "Estimate"]
  
  lm_B <- lm(y_B ~ g_B)
  beta_hat_B <- coef(summary(lm_B))["g_B", "Estimate"]
  
  y <- c(y_A, y_B)
  g <- c(g_A, g_B)
  
  lm_AUB <- lm(y ~ g)
  beta_hat_AUB <- coef(summary(lm_AUB))["g", "Estimate"]
  
  g_bar <- mean(g)
  ssg <- sum((g - g_bar) ^ 2)
  
  ssg_A <- sum((g_A - g_A_bar) ^ 2)
  ssg_B <- sum((g_B - g_B_bar) ^ 2)
  
  weight_A <- ssg_A / ssg
  weight_B <- ssg_B / ssg
  
  expected_beta_hat_AUB <- weight_A * beta_hat_A + weight_B * beta_hat_B
  
  return(
    list(
      actual = beta_hat_AUB,
      expected = expected_beta_hat_AUB
    )
  )
  
}

test <- sim_add_reg_est(
  n = 4000, m = 4000, f_A = .15, f_B = .45, beta_A = -1, beta_B = 1, sigma_A = 1, sigma_B = 1
)
