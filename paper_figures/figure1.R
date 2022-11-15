# code for figure 1
set.seed(1)
'%>%' <- magrittr::'%>%'
library(ggplot2)
source("./utils.R")
source("./get_ssg.R")

generate_fig1 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_fx_diff_noise_tradeoff_2env(
      n_e0 = 1.5e05,
      n_e1 = 1.5e05,
      se_grid = seq(from = .025, to = 1, by = .025),
      fx_diff_grid = seq(from = 0, to = 1, by = .025),
      sims_per_grid_pt = 1
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig1_df.rds")

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff_e1 = mse_e1_additive - mse_e1_GxE,
      bias_diff_e1 = bias_e1_additive ^ 2 - bias_e1_GxE ^ 2,
      var_diff_e1 = var_e1_additive - var_e1_GxE
    )

  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("MSE Additive - MSE GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  return(g1)

  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = bias_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("Bias^2 Additive - Bias^2 GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = var_diff_e1)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("Var Additive - Var GxE") +
    scale_fill_gradient2() +
    theme(plot.title = element_text(hjust = 0.5))

  gridExtra::grid.arrange(g1, g2, g3, nrow = 1)

}

generate_fig1_v2 <- function(sim_df = NULL, write_rds = FALSE) {

  if(is.null(sim_df)) {

    sim_df <- gamp::sim_var_ratio_2env(
      sigma_ratio = .884,
      se_grid = seq(from = .0001, to = .1, by = .0001),
      fx_diff_grid = seq(from = 0, to = .4, by = .0125)
    )

  }

  if(write_rds) {

    readr::write_rds(sim_df, "./rds_data/fig1_df_v2.rds")

  }

  sim_df <- sim_df %>%
    dplyr::mutate(
      mse_diff_e1 = mse_e1_additive - mse_e1_GxE,
      bias_diff_e1 = bias_e1_additive ^ 2 - bias_e1_GxE ^ 2,
      var_diff_e1 = var_e1_additive - var_e1_GxE
    )

  g1 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = mse_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("MSE Additive - MSE GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  return(g1)

  g2 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = bias_diff_e1)) +
    geom_tile() +
    scale_fill_gradient2() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("Bias^2 Additive - Bias^2 GxE") +
    theme(plot.title = element_text(hjust = 0.5))

  g3 <- ggplot(data = sim_df, aes(x = fx_diff, y = se, fill = var_diff_e1)) +
    geom_tile() +
    coord_fixed() +
    ylab("Environment Specific SE") +
    xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
    ggtitle("Var Additive - Var GxE") +
    scale_fill_gradient2() +
    theme(plot.title = element_text(hjust = 0.5))

  gridExtra::grid.arrange(g1, g2, g3, nrow = 1)

}

# first time simulation
#generate_fig1(write_rds = TRUE)
generate_fig1_v2(write_rds = T)


# subsequent analysis
sim_df <- readr::read_rds("./rds_data/fig1_df_v2.rds")

g <- generate_fig1(sim_df)

male_df <- readr::read_tsv(
  file = '../data/male_all.height.glm.linear',
  col_select = c("ID", "BETA", "SE", "OBS_CT", "P", "#CHROM", "POS")
) %>%
  dplyr::rename(
    b_m = BETA, se_m = SE, n_m = OBS_CT, p_m = P, chrom = "#CHROM", pos = POS
  )

female_df <- readr::read_tsv(
  file = '../data/female_all.height.glm.linear',
  col_select = c("ID", "BETA", "SE", "OBS_CT", "P")
) %>%
  dplyr::rename(
    b_f = BETA, se_f = SE, n_f = OBS_CT, p_f = P
  )

snp_df <- male_df %>%
  dplyr::inner_join(female_df, by = "ID")

rm(male_df)
rm(female_df)

snp_df_samp <- snp_df %>% dplyr::filter(p_m < 5e-8 | p_f < 5e-8)

# now, want to sample by LD block

ld_df <- readr::read_tsv(
  file = "../data/fourier_ls-all.bed"
) %>%
  dplyr::rename(
    chrom = chr
  ) %>%
  dplyr::mutate(
    chrom = as.numeric(gsub("chr", "", chrom, fixed = T))
  )

blocks <- nrow(ld_df)

ld_df <- ld_df %>% dplyr::mutate(ld_block = seq(1, blocks))

snp_df_samp <- sample_snps(snp_df_samp, ld_df, samps_per_block = 25)

file_names_ctr <- c(1:22)

all_c_df <- data.frame()

for (name_ctr in file_names_ctr) {

  c_df <- readr::read_tsv(
    glue::glue("/Users/ericweine/Downloads/ukb_imp_mfi/ukb_mfi_chr{name_ctr}_v3.txt"),
    col_select = c(2, 6, 7),
    col_names = FALSE
  )

  colnames(c_df) <- c("ID", "maf", "minor_allele")
  c_df <- c_df %>%
    dplyr::filter(ID %in% snp_df_samp$ID)

  all_c_df <- rbind(all_c_df, c_df)

}

snp_df_samp <- snp_df_samp %>%
  dplyr::left_join(all_c_df, by = "ID")

##### CODE TO MAKE CONTOUR PLOT ########

# hip circumference
dens <- MASS::kde2d(
  x = abs(snp_df_samp$b_f - snp_df_samp$b_m),
  y = snp_df_samp$se_m,
  h = .025,
  n = c(30, 30),
  lims = c(0, .25, 0.01, 0.05)
)

# height
dens <- MASS::kde2d(
  x = abs(snp_df_samp$b_f - snp_df_samp$b_m),
  y = snp_df_samp$se_m,
  h = .025,
  n = c(50, 50),
  lims = c(0, .1, 0.01, 0.05)
)


dens_df <- expand.grid(x = dens$x, y = dens$y)
dens_df$dens = as.vector(dens$z)
dens_df <- dens_df %>%
  dplyr::arrange(x, y)

sim_df <- sim_df %>%
  dplyr::arrange(fx_diff, se)

sim_df$dens <- dens_df$dens

sim_df <- sim_df %>%
  dplyr::mutate(
    mse_diff_e1 = mse_e1_additive - mse_e1_GxE,
    bias_diff_e1 = bias_e1_additive ^ 2 - bias_e1_GxE ^ 2,
    var_diff_e1 = var_e1_additive - var_e1_GxE
  )


sim_df <- sim_df %>%
  dplyr::filter(fx_diff < .15 & se < .15)

plot_df <- snp_df_samp %>%
  dplyr::filter(abs(b_m - b_f) < .3 & se_f < .2)

plot <- ggplot(data = plot_df, aes(y = se_f, x = abs(b_m - b_f))) +
  geom_point(alpha = .1, color = "blue") +
  geom_segment(aes(x = 0, y = 0, xend = .25, yend = 0.17), linetype = "dashed") +
  ylab("Female Coefficient SE") +
  xlab(latex2exp::TeX("$|\\hat{\\beta}_{1}^{(f)} - \\hat{\\beta}_{1}^{(m)}|$")) +
  ggtitle("Height Significant SNPs") +
  theme(plot.title = element_text(hjust = 0.5))

plot_df <- plot_df %>%
  dplyr::mutate(expected_y = abs(b_m - b_f) * .25)

ggExtra::ggMarginal(plot)

ggplot(data = dens_df, aes(x = x, y = y, fill = dens)) +
  geom_tile() +
  scale_fill_gradient2() +
  geom_segment(aes(x = .015, y = 0.01, xend = .075, yend = 0.05), linetype = "dashed") +
  ylab("Coefficient SE") +
  xlab(latex2exp::TeX("$|\\hat{\\beta}_{1}^{(f)} - \\hat{\\beta}_{1}^{(m)}|$")) +
  ggtitle("Height Significant SNPs") +
  theme(plot.title = element_text(hjust = 0.5))



  # height dec rule
  #geom_segment(aes(x = .015, y = 0.01, xend = .075, yend = 0.05), linetype = "dashed")


ggplot(data = sim_df)  +
  geom_tile(aes(x = fx_diff, y = se, fill = mse_diff_e1)) +
  scale_fill_gradient2(aes(x = fx_diff, y = se, fill = mse_diff_e1)) +
  geom_contour(aes(x = fx_diff, y = se, z = dens)) +
  coord_fixed() +
  ylab("Environment Specific SE") +
  xlab(latex2exp::TeX("$|\\beta_{1}^{(0)} - \\beta_{1}^{(1)}|$")) +
  ggtitle("MSE Additive - MSE GxE") +
  theme(plot.title = element_text(hjust = 0.5))
########################################


maf_quantiles <- quantile(snp_df_samp$maf, probs = seq(0, 1, length.out = 51))
ssg_df <- get_ssg_df(150000, maf_quantiles)
ssg_df <- setDT(ssg_df)
snp_df_samp <- setDT(snp_df_samp)

snp_df_samp <- snp_df_samp[ssg_df, on = .(maf >= maf_low, maf < maf_high),
              nomatch = 0, ]

snp_df_samp <- snp_df_samp %>%
  dplyr::mutate(
    sigma_m = gamp:::get_theo_obs_noise(se_m, n_m, maf),
    sigma_f = gamp:::get_theo_obs_noise(se_f, n_f, maf)
  )

get_decision_rule <- function(ssg, b_m, b_f, sigma_m, sigma_f) {

  mse_GxE <- (sigma_m ^ 2) / ssg

  bias_additive <- (b_m * ssg + b_f * ssg) / (2 * ssg) - b_m

  var_additive <- ((sigma_m ^ 2) * ssg + (sigma_f ^ 2) * ssg) / (4 * (ssg ^ 2))

  mse_additive <- bias_additive ^ 2 + var_additive

  if (mse_GxE < mse_additive) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}

# now, I want to run ash for each individual maf bin
# and then perhaps do some thresholding
# I also have to get the appropriate decision rule

under_thresh_vec <- c()

thresh_vec <- seq(from = 0, to = .5, by = .025)

for (i in c(2:length(thresh_vec))) {

  print(thresh_vec[i])

  thresh_df <- snp_df_samp %>%
    dplyr::filter(maf <= thresh_vec[i] & maf > thresh_vec[i - 1])

  # print(glue::glue("Running ash with {nrow(thresh_df)} snps"))
  #
  # male_ash <- ashr::ash(
  #   betahat = thresh_df$b_m,
  #   sebetahat = thresh_df$se_m,
  #   mixcompdist = "normal",
  #   outputlevel = 2
  # )
  # female_ash <- ashr::ash(
  #   betahat = thresh_df$b_f,
  #   sebetahat = thresh_df$se_f,
  #   mixcompdist = "normal",
  #   outputlevel = 2
  # )
  #
  # thresh_df <- thresh_df %>%
  #   dplyr::mutate(
  #     b_m_pm = ashr::get_pm(male_ash),
  #     b_m_lfsr = ashr::get_lfsr(male_ash),
  #     b_f_pm = ashr::get_pm(female_ash),
  #     b_f_lfsr = ashr::get_lfsr(female_ash)
  #   )
  #
  # thresh_df <- thresh_df %>%
  #   dplyr::filter(b_m_lfsr < .05 | b_f_lfsr < .05)

  under_thresh_ct <- 0

  for(j in 1:nrow(thresh_df)) {

    under_thresh_ct <- under_thresh_ct + as.numeric(
      get_decision_rule(
        thresh_df$ssg[i], thresh_df$b_m[i],thresh_df$b_f[i], 3.682341, 0.6177109
      )
    )

  }

  #print(glue::glue("{nrow(thresh_df)} remaining after lfsr filtering"))

  # now, apply the decision rule here

  pct_under_thresh <- under_thresh_ct / nrow(thresh_df)
  under_thresh_vec <- c(under_thresh_vec, pct_under_thresh)

}


# male_ash <- ashr::ash(
#   betahat = snp_df_samp$b_m,
#   sebetahat = snp_df_samp$se_m,
#   mixcompdist = "normal",
#   outputlevel = 1
#   )
# female_ash <- ashr::ash(
#   betahat = snp_df_samp$b_f,
#   sebetahat = snp_df_samp$se_f,
#   mixcompdist = "normal",
#   outputlevel = 1
#   )
#
# snp_df_samp <- snp_df_samp %>%
#   dplyr::mutate(
#     b_m_pm = ashr::get_pm(male_ash),
#     b_f_pm = ashr::get_pm(female_ash),
#   )
#
# # here, I can make a plot where I subset based on MAF
# # and then look at the percentage of snps in that bucket
# # that goes over the decision boundary
#
#
#
# #snp_df_samp <- snp_df_samp %>%
# #  dplyr::filter(se_m <= 1)
#
# snp_df_samp <- snp_df_samp %>%
#   dplyr::mutate(fx_diff = abs(b_m_pm - b_f_pm))
#
# #g2 <- ggplot2::ggplot(data = snp_df_samp, aes(x = fx_diff, y = se_m)) +
# #  ggplot2::geom_point() +
# #  ggplot2::coord_fixed() +
# #  ggplot2::xlim(0, 1) +
# #  ggplot2::ylim(0, 1)
#
# #gridExtra::grid.arrange(g, ggExtra::ggMarginal(g2), nrow = 1)
#
# under_thresh_vec <- c()
#
# thresh_vec <- seq(from = 0, to = .5, by = .025)
#
# for (i in c(2:length(thresh_vec))) {
#
#   thresh_df <- snp_df_samp %>%
#     dplyr::filter(maf <= thresh_vec[i] & maf > thresh_vec[i - 1])
#
#   num_under_thresh <- sum(thresh_df$se_m < thresh_df$expected_se)
#
#   pct_under_thresh <- num_under_thresh / nrow(thresh_df)
#   under_thresh_vec <- c(under_thresh_vec, pct_under_thresh)
#
# }

out_df <- data.frame(
  maf_thresh = seq(from = 0.025, to = .5, by = .025),
  pct_under = under_thresh_vec
)

ggplot(data = out_df, aes(x = maf_thresh, y = pct_under)) +
  geom_point() +
  ylab("Proportion Below Decision Boundary") +
  xlab("MAF Bucket") +
  ggtitle("Hip Circ")


maf_out <- function(y) {

  sol1 <- (-y^2 + sqrt(y^4 - 4 * y^2)) / (-2 * y^2)
  sol2 <- (-y^2 - sqrt(y^4 - 4 * y^2)) / (-2 * y^2)
  return(c(sol1, sol2))

}
