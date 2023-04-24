# Here, I want to roughly simulate a process similar to Traglia
# then I can see the results and demonstrate why their method is bad

Sigma <- list()

Sigma[[1]] <- .0025 * gamp:::make_amp_cov_mat(
  desired_corr = 1, amp_coef = 3, amp_e0 = FALSE
)

simulated_fx <- gamp:::generate_mash_fx(
  n = 9600000, p = 2, Sigma = Sigma, pi = c(1), s = 0
)


# I can do LD pruning here later as well
library(dplyr)

male_df <- readr::read_tsv(
  file = '/Users/ericweine/Downloads/diastolicBP_auto/male_all.diastolicBP_auto.glm.linear',
  col_select = c("ID", "BETA", "SE", "P")
) %>%
  dplyr::rename(
    b_m = BETA, se_m = SE, p_m = P
  )

female_df <- readr::read_tsv(
  file = '/Users/ericweine/Downloads/diastolicBP_auto/female_all.diastolicBP_auto.glm.linear',
  col_select = c("ID", "BETA", "SE", "P")
) %>%
  dplyr::rename(
    b_f = BETA, se_f = SE, p_f = P
  )

snp_df <- male_df %>%
  dplyr::inner_join(female_df, by = "ID")

rm("male_df")
rm("female_df")
snp_df <- snp_df %>%
  dplyr::distinct(ID, .keep_all = TRUE)

# snp_df <- snp_df %>%
#   dplyr::filter(
#     p_m < 5e-8 | p_f < 5e-8
#   )

# snp_df <- snp_df %>%
#   dplyr::filter(
#     se_m < 3 & se_f < 3
#   )

# set.seed(1)
# snp_samp_df <- snp_df %>%
#   dplyr::sample_n(10000)

# library(ggplot2)
# 
# ggplot(snp_df, aes(x = b_m, y = b_f)) +
#   geom_point(alpha = .2)

# now, I should select a random subset

# I should also conduct this test on the real data and examine the results
# it's hard to judge because mash is good at filtering noise out
# which these methods may be worse at
# however, that's another way that maybe I can get the sim to look better
# I'm not totally sure

snp_df <- snp_df %>%
  dplyr::mutate(
    sex_t_stat = abs(b_m - b_f) / sqrt((se_m ^ 2) + (se_f ^ 2))
  )

snp_df <- snp_df %>%
  dplyr::mutate(
    sex_diff_p_val = 2 * (1 - pnorm(sex_t_stat))
  )

snp_sex_diff_sig <- snp_df %>%
  dplyr::filter(sex_diff_p_val < 5e-8)

snp_df <- snp_df %>%
  select(c(se_f, se_m))

#e0_sds <- rchisq(nrow(simulated_fx), 1, 1)
sd_idx <- sample.int(n = nrow(snp_df), size = 9600000)

e0_sds <- snp_df$se_m[sd_idx]

e0_sim_mles <- rnorm(
  n = nrow(simulated_fx), 
  mean = simulated_fx[, 1], 
  sd = e0_sds
)

e0_sd_est <- sqrt(((e0_sds ^ 2) / (150000)) * rchisq(9600000, 150000))

e1_sds <- snp_df$se_f[sd_idx]

e1_sim_mles <- rnorm(
  n = nrow(simulated_fx), 
  mean = simulated_fx[, 2], 
  sd = e1_sds
)

e1_sd_est <- sqrt(((e1_sds ^ 2) / (150000)) * rchisq(9600000, 150000))

prior_sim_df <- data.frame(
  male_fx = simulated_fx[, 1],
  female_fx = simulated_fx[, 2],
  male_mle = e0_sim_mles,
  female_mle = e1_sim_mles,
  male_sd = e0_sds,
  female_sd = e1_sds,
  male_sd_est = e0_sd_est,
  female_sd_est = e1_sd_est,
  mix_comp = simulated_fx[, 3]
)

prior_sim_df <- prior_sim_df %>%
  dplyr::mutate(
    sex_diff_stat = abs(male_mle - female_mle) / sqrt((male_sd_est ^ 2) + (female_sd_est ^ 2))
  )

prior_sim_df <- prior_sim_df %>%
  dplyr::mutate(
    sex_diff_p = 2 * (1 - pnorm(sex_diff_stat))
  )

sig_df <- prior_sim_df %>%
  dplyr::filter(sex_diff_p < 5e-8)

sig_df2 <- prior_sim_df %>%
  dplyr::filter(sex_diff_p < 5e-5)

# now, I should be able to run analysis on this and see what the results look like
# I'm not sure how much this will really look like the data in question
# I will probably have to ask Arbel about that

# first, plot the simulated data
library(ggplot2)
ggplot(data = prior_sim_df, aes(x = male_mle, y = female_mle)) +
  geom_point(alpha = .2)

ggplot(data = sig_df, aes(x = male_mle, y = female_mle)) +
  geom_point(alpha = .2)
