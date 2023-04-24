set.seed(1)
data_dir <- "/Users/ericweine/Downloads/doi_10.5061_dryad.hhmgqnkgr__v7"
library(dplyr)
summary_table <- read.delim(
  glue::glue('{data_dir}/SummaryTable_allsites_12Nov20.txt')
)

sites_df <- data.frame(stringr::str_split_fixed(summary_table$site, ":", 2))
colnames(sites_df) <- c("chromosome", "site_id")
sites_df <- sites_df %>%
  dplyr::mutate(site_id = as.numeric(site_id))
# split into blocks of a certain length
split_into_LD_blocks <- function(df, block_length) {
  block_range <- seq(from = min(df$site_id), to = max(df$site_id), by = block_length)
  df %>%
    dplyr::mutate(block_id = plyr::laply(site_id, function(x) sum(x > block_range)))
}
# group by chromosome and then split into blocks
sites_df <- sites_df %>%
  dplyr::group_by(chromosome) %>%
  dplyr::group_modify(~ split_into_LD_blocks(.x, 1e4))
# Now, want to sample one SNP from each group
sites_sample_df <- sites_df %>%
  dplyr::ungroup() %>%
  dplyr::sample_frac() %>% #randomly shuffle df
  dplyr::distinct(chromosome, block_id, .keep_all = TRUE) %>%
  dplyr::select(chromosome, site_id)
# Reconstruct site names
selected_sites <- purrr::pmap_chr(
  list(sites_sample_df$chromosome, sites_sample_df$site_id),
  function(x, y) glue::glue("{x}:{y}")
)
summary_table_samp <- summary_table %>%
  dplyr::filter(site %in% selected_sites) %>%
  dplyr::select(c(site, pval_CTRL, pval_HS, coef_CTRL, coef_HS, sig_cat))
# replace 0 p-values with small numbers
summary_table_samp <- summary_table_samp %>%
  dplyr::mutate(
    pval_CTRL = pmax(1e-8, pval_CTRL),
    pval_HS = pmax(1e-8, pval_HS)
  )
# construct std error estimates from coefficients and p-values
summary_table_samp <- summary_table_samp %>%
  dplyr::mutate(
    std_error_ctrl = abs(coef_CTRL) / qnorm((2 - pval_CTRL) / 2),
    std_error_hs = abs(coef_HS) / qnorm((2 - pval_HS) / 2)
  )

reg_fx_mat <- t(matrix(
  data = c(summary_table_samp$coef_CTRL, summary_table_samp$coef_HS),
  nrow = 2,
  byrow = TRUE
))
colnames(reg_fx_mat) <- c("ctrl", "hs")
reg_se_mat <- t(matrix(
  data = c(summary_table_samp$std_error_ctrl, summary_table_samp$std_error_hs),
  nrow = 2,
  byrow = TRUE
))
colnames(reg_se_mat) <- c("ctrl", "hs")
mash_data <- mashr::mash_set_data(reg_fx_mat, reg_se_mat)
# Now, want to construct covariance matrices to feed into mash
cov_mat_list <- list()
cov_mat_list[['hs_spec']] <- matrix(
  data = c(0, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(
    rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
  )
)
cov_mat_list[['ctrl_spec']] <- matrix(
  data = c(1, 0, 0, 0), nrow = 2, byrow = TRUE, dimnames = list(
    rows = c("ctrl", "hs"), cols = c("ctrl", "hs")
  )
)
desired_corrs <- seq(from = -1, to = 1, by = (1/3))
desired_amp <- c(2, 1.5, 1.25)
for(corr in desired_corrs) {
  cov_mat_list[[glue::glue('equal_corr_{corr}')]] <- gamp:::make_amp_cov_mat(
    desired_corr = corr, amp = FALSE
  )
  for(cond in c("hs", "ctrl")) {
    for(amp in desired_amp) {
      cov_mat_list[[glue::glue('{cond}_amp_{amp}_corr_{corr}')]] <- gamp:::make_amp_cov_mat(
        desired_corr = corr, amp_e0 = (cond == "ctrl"), amp_coef = amp
      )
    }
  }
}
mash_out <- mashr::mash(
  data = mash_data,
  Ulist = cov_mat_list,
  algorithm.version = "Rcpp"
)


# here, I have to make sure that I get all of the data

mash_out_all <- mashr::mash(
  data = mash_data,
  fixg = TRUE,
  g = ashr::get_fitted_g(mash_out),
  algorithm.version = "Rcpp",
  outputlevel = 3
)

cov_mat_ests <- mashr::get_estimated_pi(mash_out)

# mash_post_samp <- mashr::mash(
#   data = mash_data,
#   algorithm.version = "R",
#   posterior_samples = 1,
#   fixg = TRUE,
#   g = ashr::get_fitted_g(mash_out)
# )
# 
# samp <- mashr::get_samples(mash_post_samp)[,,"sample_1"]
# 
# # Now, I can use a posterior sample from this
# # and add some noise
# # to generate one sample of what the MLEs of the data look like
# # then I can run the Pallares algorithm
# # and show all of the discoveries you would make. 
# 
# sim_data_ctrl <- rnorm(
#   n = nrow(samp), 
#   mean = samp[,"ctrl"], 
#   sd = reg_se_mat[, "ctrl"]
# )
# 
# sim_data_hs <- rnorm(
#   n = nrow(samp), 
#   mean = samp[,"hs"], 
#   sd = reg_se_mat[, "hs"]
# )
# 
# sim_df <- data.frame(
#   ctrl = sim_data_ctrl,
#   hs = sim_data_hs
# )
# 
# ggplot(data = sim_df, aes(x = ctrl, y = hs)) +
#   geom_point()
# 
# ggplot(data = summary_table_samp, aes(x = coef_CTRL, y = coef_HS)) +
#   geom_point()
# 
# # I suppose that it actually probably makes more sense
# # to simulate effects from the fitted mash prior
# # then, I know which effect was which
# 
# # I will prepare a relatively simple mixture distribution
# # then, I can run the pallares procedure on it

fitted_g <- ashr::get_fitted_g(mash_out)
Sigma <- list()
#Sigma[[1]] <- fitted_g$Ulist$equal_corr_1 * .04
Sigma[[1]] <- .01 * gamp:::make_amp_cov_mat(
  desired_corr = 1, amp_coef = 1.4, amp_e0 = FALSE
)

mu <- list()
mu[[1]] <- c(-.125, -.15)

#Sig_pi <- c(.55, .35) / sum(c(.55, .35))
Sig_pi <- 1

simulated_fx <- gamp:::generate_mash_fx(
  n = 50000, p = 2, Sigma = Sigma, pi = Sig_pi, s = .6, mu = mu
)

sd_idx <- sample.int(nrow(summary_table_samp), size = 50000, replace = TRUE)

#e0_sds <- invgamma::rinvgamma(nrow(simulated_fx), shape = 6, scale = 1.5)
# e0_sds <- sample(
#   x = summary_table_samp$std_error_ctrl,
#   replace = TRUE,
#   size = 50000
# )

e0_sds <- summary_table_samp$std_error_ctrl[sd_idx]

e0_sim_mles <- rnorm(
  n = nrow(simulated_fx), 
  mean = simulated_fx[, 1], 
  sd = e0_sds
)

#e1_sds <- invgamma::rinvgamma(nrow(simulated_fx), shape = 6, scale = 1.5)
# e1_sds <- sample(
#   x = summary_table_samp$std_error_hs,
#   replace = TRUE,
#   size = 50000
# )

e1_sds <- summary_table_samp$std_error_hs[sd_idx]

e1_sim_mles <- rnorm(
  n = nrow(simulated_fx), 
  mean = simulated_fx[, 2], 
  sd = e1_sds
)

# now, put everything together in a matrix
# then, I can run an approximation to the analysis

prior_sim_df <- data.frame(
  e0_fx = simulated_fx[, 1],
  e1_fx = simulated_fx[, 2],
  e0_mle = e0_sim_mles,
  e1_mle = e1_sim_mles,
  e0_sd = e0_sds,
  e1_sd = e1_sds,
  mix_comp = simulated_fx[, 3]
)

# ggplot(data = prior_sim_df, aes(x = e0_mle, y = e1_mle)) +
#   geom_point()

# this data looks reasonable
# Now, the goal is to run the same analysis as Pallares

prior_sim_df <- prior_sim_df %>%
  dplyr::mutate(
    e0_pval = 2*pnorm(-abs(e0_mle / e0_sd)),
    e1_pval = 2*pnorm(-abs(e1_mle / e1_sd)),
  )

prior_sim_df <- prior_sim_df %>%
  dplyr::mutate(
    e0_qval = qvalue::qvalue(e0_pval)$qvalues,
    e1_qval = qvalue::qvalue(e1_pval)$qvalues
  )

# Now, I must classify the signals
prior_sim_df <- prior_sim_df %>%
  dplyr::mutate(
    sig_class = dplyr::case_when(
      e0_qval > .01 & e1_qval > .01 ~ "null",
      (e0_qval < .01 & e1_pval < .1) | (e1_qval < .1 & e0_pval < .1) ~ "shared",
      e0_qval < .01 & e1_pval >= .1 ~ "CTRL",
      e1_qval < .01 & e0_pval >= .1 ~ "HS",
      TRUE ~ "error"
    )
  )

signif_df <- prior_sim_df %>%
  dplyr::filter(sig_class != "null")

prior_sim_samp_df <- prior_sim_df %>%
  dplyr::sample_n(12030)

library(ggplot2)

sim_data_plot <- ggplot(data = prior_sim_samp_df, aes(x = e0_mle, y = e1_mle)) +
  geom_point(alpha = .1) +
  xlab("Effect Estimate in Control Group") +
  ylab("Effect Estimate in High Sugar Group") +
  ggtitle("Simulated Data from\nPervasive Amplification Model") +
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    plot.tag.position=c(1.25, 0.5)
  ) +
  scale_x_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  labs(tag = "Random Variants")

sim_sig_plot <- ggplot(data = signif_df, aes(x = e0_mle, y = e1_mle, color = sig_class)) +
  geom_point(alpha = .175) +
  xlab("Effect Estimate in Control Group") +
  ylab("Effect Estimate in High Sugar Group") +
  guides(color=guide_legend(title="Genetic Effect Classification", override.aes = list(alpha = .5))) +
  scale_color_hue(
    labels = c(
      "Specific to Control Diet", 
      "Specific to High Sugar Diet", 
      "Shared between Diets"
    )
  ) +
  scale_x_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.box.margin=margin(l=20),
    plot.tag.position=c(1.25, 0.5),
    axis.title.y.right = element_text(size=5)
  ) +
  labs(tag = "Ascertained as Significant")

real_data_plot <- ggplot(data = summary_table_samp, aes(x = coef_CTRL, y = coef_HS)) +
  geom_point(alpha = .1) +
  xlab("Effect Estimate in Control Group") +
  ylab("Effect Estimate in High Sugar Group") +
  ggtitle("Genetic Effect on\nLongevity in Flies (Pallares et al., 2023)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.box.margin=margin(l=20)
  ) +
  scale_x_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits=c(-1, .75), expand = c(0.01, 0.01))

summary_table_sig <- summary_table %>%
  dplyr::filter(sig_cat != "NS") %>%
  dplyr::filter(coef_CTRL < 1 & coef_HS < 1)

real_sig_plot <- ggplot(data = summary_table_sig, aes(x = coef_CTRL, y = coef_HS, color = sig_cat)) +
  geom_point(alpha = .175) +
  xlab("Effect Estimate in Control Group") +
  ylab("Effect Estimate in High Sugar Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8)
  ) +
  guides(
    color=guide_legend(title="Genetic Effect Classification", override.aes = list(alpha = .5))) +
  scale_color_hue(
    labels = c(
      "Specific to Control Diet", 
      "Specific to High Sugar Diet", 
      "Shared between Diets"
    )
  ) +
  scale_x_continuous(limits=c(-1, .75), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits=c(-1, .75), expand = c(0.01, 0.01))
  

library(ggpubr)

fig4 <- ggpubr::ggarrange(
  real_data_plot + rremove("ylab") + rremove("xlab"), 
  sim_data_plot + rremove("ylab") + rremove("xlab"), 
  real_sig_plot + rremove("ylab") + rremove("xlab"), 
  sim_sig_plot + rremove("xlab") + rremove("ylab"),
  common.legend = T,
  legend = "right",
  labels = "AUTO"
)

library(grid)

annotate_figure(fig4, 
                left = textGrob("Effect Estimate in Flies Given High Sugar Diet", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 11)),
                bottom = textGrob("Effect Estimate in Flies Given Control Diet", gp = gpar(cex = 1, fontsize = 11), hjust = .5)
                )

