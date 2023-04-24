# first, I want to read in the data from lonestar
# then I can make a grid
'%>%' <- magrittr::'%>%'
set.seed(1)
make_dec_rule_male_df <- function(
  data_dir, trait, n, m, n_samp, pval_thresh = 1
) {
  
  cv_df <- readr::read_csv(
    file = glue::glue("{data_dir}/cv_{trait}_fold1_ash_full.csv")
  )
  
  cv_df <- cv_df %>%
    dplyr::mutate(p_val = 2 * pnorm(
      cv_df$b_train_male_1 / cv_df$se_train_male_1, 
      lower.tail = FALSE
      )
    )
  
  cv_df <- cv_df %>%
    dplyr::filter(p_val < pval_thresh)
  
  dec_rule <- gamp:::get_sigma_ratio_slope(
    median(cv_df$se_train_male_1 / cv_df$se_train_female_1), n = n, m = m
  )
  
  male_cv_df <- cv_df %>%
    dplyr::select(se_train_male_1, expected_add_male_bias, actual_male_add_mse_advantage, b_train_male_1) %>%
    dplyr::mutate(fx_diff = 2 * abs(expected_add_male_bias)) %>%
    dplyr::rename(se = se_train_male_1, mse_diff = actual_male_add_mse_advantage) %>%
    dplyr::select(-expected_add_male_bias) %>%
    dplyr::mutate(expected_se = dec_rule * fx_diff)
  
  pct_under_line <- nrow(male_cv_df %>% dplyr::filter(se < expected_se)) / nrow(male_cv_df)
  
  if (nrow(male_cv_df) > n_samp) {
    
    samp_cv_df <- male_cv_df %>%
      dplyr::sample_n(n_samp)
    
  } else {
    
    samp_cv_df <- male_cv_df
    
  }
  
  return(
    list(
      slope = dec_rule,
      cv_df = samp_cv_df,
      pct_under_line = pct_under_line
    )
  )
  
}

make_dec_rule_female_df <- function(
    data_dir, trait, n, m, n_samp, pval_thresh = 1
) {
  
  cv_df <- readr::read_csv(
    file = glue::glue("{data_dir}/cv_{trait}_fold1_ash_full.csv")
  )
  
  cv_df <- cv_df %>%
    dplyr::mutate(p_val = 2 * pnorm(
      cv_df$b_train_female_1 / cv_df$se_train_female_1, 
      lower.tail = FALSE
    )
    )
  
  cv_df <- cv_df %>%
    dplyr::filter(p_val < pval_thresh)
  
  dec_rule <- gamp:::get_sigma_ratio_slope(
    median(cv_df$se_train_female_1 / cv_df$se_train_fefemale_1), n = n, m = m
  )
  
  female_cv_df <- cv_df %>%
    dplyr::select(se_train_female_1, expected_add_female_bias, actual_female_add_mse_advantage, b_train_female_1) %>%
    dplyr::mutate(fx_diff = 2 * abs(expected_add_female_bias)) %>%
    dplyr::rename(se = se_train_female_1, mse_diff = actual_female_add_mse_advantage) %>%
    dplyr::select(-expected_add_female_bias) %>%
    dplyr::mutate(expected_se = dec_rule * fx_diff)
  
  pct_under_line <- nrow(female_cv_df %>% dplyr::filter(se < expected_se)) / nrow(female_cv_df)
  
  if (nrow(female_cv_df) > n_samp) {
    
    samp_cv_df <- female_cv_df %>%
      dplyr::sample_n(n_samp)
    
  } else {
    
    samp_cv_df <- female_cv_df
    
  }
  
  return(
    list(
      slope = dec_rule,
      cv_df = samp_cv_df,
      pct_under_line = pct_under_line
    )
  )
  
}

make_crossval_plot_male <- function(data_dir, trait) {

  cv_df <- readr::read_csv(
    file = glue::glue("{data_dir}/cv_{trait}_fold1_ash_full.csv")
  )

  male_cv_df <- cv_df %>%
    dplyr::filter(
      expected_male_add_mse_advantage > quantile(expected_male_add_mse_advantage, probs = .01) &
        expected_male_add_mse_advantage < quantile(expected_male_add_mse_advantage, probs = .99)
    )
  
  male_cv_df <- male_cv_df %>%
    dplyr::mutate(
      percentile = cut(
        expected_male_add_mse_advantage,
        seq(min(expected_male_add_mse_advantage), max(expected_male_add_mse_advantage), length.out = 10),
        include.lowest = TRUE,
        labels = FALSE
      )
    )
  
  male_summ_df <- male_cv_df %>%
    dplyr::group_by(percentile) %>%
    dplyr::summarize(
      weighted_dec_advantage_add_male = weighted.mean(actual_male_add_mse_advantage, 1 / (se_test_male_1 ^ 2)),
      unweighted_dec_advantage_add_male = mean(actual_male_add_mse_advantage),
      weighted_sd_advantage_add_male = sqrt(Hmisc::wtd.var(actual_male_add_mse_advantage, 1 / (se_test_male_1 ^ 2))),
      expected_advantage = mean(expected_male_add_mse_advantage),
      samp_size = dplyr::n(),
      sd = var(actual_male_add_mse_advantage)
    )
  
  return(male_summ_df)

}

make_crossval_plot_female <- function(data_dir, trait) {
  
  cv_df <- readr::read_csv(
    file = glue::glue("{data_dir}/cv_{trait}_fold1_ash_full.csv")
  )
  
  female_cv_df <- cv_df %>%
    dplyr::filter(
      expected_female_add_mse_advantage > quantile(expected_female_add_mse_advantage, probs = .01) &
        expected_female_add_mse_advantage < quantile(expected_female_add_mse_advantage, probs = .99)
    )
  
  female_cv_df <- female_cv_df %>%
    dplyr::mutate(
      percentile = cut(
        expected_female_add_mse_advantage,
        seq(min(expected_female_add_mse_advantage), max(expected_female_add_mse_advantage), length.out = 10),
        include.lowest = TRUE,
        labels = FALSE
      )
    )
  
  female_summ_df <- female_cv_df %>%
    dplyr::group_by(percentile) %>%
    dplyr::summarize(
      weighted_dec_advantage_add_female = weighted.mean(actual_female_add_mse_advantage, 1 / (se_test_female_1 ^ 2)),
      unweighted_dec_advantage_add_female = mean(actual_female_add_mse_advantage),
      weighted_sd_advantage_add_female = sqrt(Hmisc::wtd.var(actual_female_add_mse_advantage, 1 / (se_test_female_1 ^ 2))),
      expected_advantage = mean(expected_female_add_mse_advantage),
      samp_size = dplyr::n(),
      sd = var(actual_female_add_mse_advantage)
    )
  
  return(female_summ_df)
  
}

library(ggplot2)

igf1_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "IGF1",
  n = 148128,
  m = 171544,
  n_samp = 100000
)

g_igf1 <- ggplot(data = igf1_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .2, yend = .2 * igf1_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (nmol/L)") +
  ylab("Standard Error (nmol/L)") +
  ggtitle("IGF1 Random Variants") +
  xlim(0, .2) +
  ylim(0, 4)

igf1_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "IGF1",
  n = 148128,
  m = 171544,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_igf1_top_hits <- ggplot(data = igf1_info_top_hits$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .25, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .2, yend = .2 * igf1_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (nmol/L)") +
  ylab("Standard Error (nmol/L)") +
  ggtitle("IGF1 Top Hits (p < 5e-8 in Males)") +
  xlim(0, .2) +
  ylim(0, (1/3)) 

library(ggpubr)

png(
  file = "~/Documents/paper_gwas_amplification/images/fig3.png",
  width = 4000,
  height = 2000,
  res = 300
)

igf1_fig <- ggarrange(
  g_igf1 + rremove("ylab") + rremove("xlab"),
  g_igf1_top_hits + rremove("ylab") + rremove("xlab"),
  labels = "AUTO"
)

annotate_figure(
  igf1_fig,
  left = textGrob("Standard Error (nmol/L)", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 14)),
  bottom = textGrob("Difference in Context Specific Effects (nmol/L)", gp = gpar(cex = 1, fontsize = 14), hjust = .575))
dev.off()

albumin_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "albumin",
  n = 137529,
  m = 156809,
  n_samp = 100000
)

albumin_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "albumin",
  n = 137529,
  m = 156809,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_albumin <- ggplot(data = albumin_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .02, yend = .02 * albumin_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (g/L)") +
  ylab("Standard Error (g/L)") +
  ggtitle("Albumin") +
  xlim(0, .02) +
  ylim(0, 1.5)

fvc_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "FVC_best",
  n = 119821,
  m = 135605,
  n_samp = 100000
)

fvc_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "FVC_best",
  n = 119821,
  m = 135605,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_fvc <- ggplot(data = fvc_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .004, yend = .004 * fvc_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (L)") +
  ylab("Standard Error (L)") +
  ggtitle("FVC") +
  xlim(0, .004) +
  ylim(0, .75)

hba1c_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "HbA1c",
  n = 148771,
  m = 172497,
  n_samp = 100000
)

hba1c_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "HbA1c",
  n = 148771,
  m = 172497,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_hba1c <- ggplot(data = hba1c_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .015, yend = .015 * hba1c_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmol/mol)") +
  ylab("Standard Error (mmol/mol)") +
  ggtitle("HbA1C") +
  xlim(0, .015) +
  ylim(0, 5)

rbc_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "RBC_count",
  n = 151792,
  m = 175344,
  n_samp = 100000
)

rbc_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "RBC_count",
  n = 151792,
  m = 175344,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_rbc <- ggplot(data = rbc_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .0025, yend = .0025 * rbc_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (10^12 cells/L)") +
  ylab("Standard Error (10^12 cells/L)") +
  ggtitle("RBC") +
  xlim(0, .0025) +
  ylim(0, .5)

shbg_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "SHBG",
  n = 136468,
  m = 155040,
  n_samp = 100000
)

shbg_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "SHBG",
  n = 136468,
  m = 155040,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_shbg <- ggplot(data = shbg_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 1.25, yend = 1.25 * shbg_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (nmol/L)") +
  ylab("Standard Error (nmol/L)") +
  ggtitle("SHBG") +
  xlim(0, 1.25) +
  ylim(0, 15)

arm_fatfree_mass_L_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "arm_fatfree_mass_L",
  n = 153102,
  m = 177972,
  n_samp = 100000
)

arm_fatfree_mass_L_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "arm_fatfree_mass_L",
  n = 153102,
  m = 177972,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_arm_L <- ggplot(data = arm_fatfree_mass_L_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .025, yend = .025 * arm_fatfree_mass_L_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (Kg)") +
  ylab("Standard Error (Kg)") +
  ggtitle("Left Arm Fat-Free Mass") +
  xlim(0, .025) +
  ylim(0, .6)

arm_fatfree_mass_R_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "arm_fatfree_mass_R",
  n = 153135,
  m = 178001,
  n_samp = 100000
)

arm_fatfree_mass_R_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "arm_fatfree_mass_R",
  n = 153135,
  m = 178001,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_arm_R <- ggplot(data = arm_fatfree_mass_R_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .025, yend = .025 * arm_fatfree_mass_R_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (Kg)") +
  ylab("Standard Error (Kg)") +
  ggtitle("Right Arm Fat-Free Mass") +
  xlim(0, .025) +
  ylim(0, .6)

bmi_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "bmi",
  n = 155556,
  m = 180464,
  n_samp = 100000
)

bmi_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "bmi",
  n = 155556,
  m = 180464,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_bmi <- ggplot(data = bmi_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .05, yend = .05 * bmi_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (Kg/m^2)") +
  ylab("Standard Error (Kg/m^2)") +
  ggtitle("BMI") +
  xlim(0, .05) +
  ylim(0, 3)

calcium_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "calcium",
  n = 137489,
  m = 156741,
  n_samp = 100000
)

calcium_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "calcium",
  n = 137489,
  m = 156741,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_calcium <- ggplot(data = calcium_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .0005, yend = .0005 * calcium_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmol/L)") +
  ylab("Standard Error (mmol/L)") +
  ggtitle("Calcium") +
  xlim(0, .0005) +
  ylim(0, .1)

creatinine_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "creatinine",
  n = 148824,
  m = 172409,
  n_samp = 100000
)

creatinine_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "creatinine",
  n = 148824,
  m = 172409,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_creatinine <- ggplot(data = creatinine_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 5, yend = 5 * creatinine_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (umol/L)") +
  ylab("Standard Error (umol/L)") +
  ggtitle("Creatinine") +
  ylim(0, 20) +
  xlim(0, 5)

diastolicBP_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "diastolicBP_auto",
  n = 145992,
  m = 168906,
  n_samp = 100000
)

diastolicBP_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "diastolicBP_auto",
  n = 145992,
  m = 168906,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_diastolicBP <- ggplot(data = diastolicBP_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .1, yend = .1 * diastolicBP_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmHg)") +
  ylab("Standard Error (mmHg)") +
  ggtitle("Diastolic BP")

eosinophil_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "eosinophil_perc",
  n = 151524,
  m = 175049,
  n_samp = 100000
)

eosinophil_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "eosinophil_perc",
  n = 151524,
  m = 175049,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_eosinophil <- ggplot(data = eosinophil_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .0125, yend = .0125 * eosinophil_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (cells/L)") +
  ylab("Standard Error (cells/L)") +
  ggtitle("Eosinophil") +
  xlim(0, .0125) +
  ylim(0, 1.75)

height_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "height",
  n = 155732,
  m = 180655,
  n_samp = 100000
)

height_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "height",
  n = 155732,
  m = 180655,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_height <- ggplot(data = height_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .1, yend = .1 * height_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (cm)") +
  ylab("Standard Error (cm)") +
  ggtitle("Height") +
  xlim(0, .1) + 
  ylim(0, 5)

hip_circ_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "hip_circ",
  n = 155828,
  m = 180685,
  n_samp = 100000
)

hip_circ_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "hip_circ",
  n = 155828,
  m = 180685,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_hip_circ <- ggplot(data = hip_circ_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .22, yend = .22 * hip_circ_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (cm)") +
  ylab("Standard Error (cm)") +
  ggtitle("Hip Circumference") +
  ylim(0, 5) +
  xlim(0, .22)

lymphocyte_perc_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "lymphocyte_perc",
  n = 151524,
  m = 175049,
  n_samp = 100000
)

lymphocyte_perc_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "lymphocyte_perc",
  n = 151524,
  m = 175049,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_lymphocyte_perc <- ggplot(data = lymphocyte_perc_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .022, yend = .022 * lymphocyte_perc_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (%)") +
  ylab("Standard Error (%)") +
  ggtitle("Lymphocyte %") +
  xlim(0, .022) +
  ylim(0, 7)

protein_total_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "protein_total",
  n = 137360,
  m = 156652,
  n_samp = 100000
)

protein_total_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "protein_total",
  n = 137360,
  m = 156652,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_protein_total <- ggplot(data = protein_total_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .025, yend = .025 * protein_total_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmol/L)") +
  ylab("Standard Error (mmol/L)") +
  ggtitle("Protein Total") +
  xlim(0, .025) +
  ylim(0, 4.5)

pulse_rate_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "pulse_rate",
  n = 145992,
  m = 168906,
  n_samp = 100000
)

pulse_rate_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "pulse_rate",
  n = 145992,
  m = 168906,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_pulse_rate <- ggplot(data = pulse_rate_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .05, yend = .05 * pulse_rate_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (bpm)") +
  ylab("Standard Error (bpm)") +
  ggtitle("Pulse Rate") +
  xlim(0, .05) +
  ylim(0, 7.5)

systolicBP_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "systolicBP_auto",
  n = 145989,
  m = 168903,
  n_samp = 100000
)

systolicBP_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "systolicBP_auto",
  n = 145989,
  m = 168903,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_systolicBP <- ggplot(data = systolicBP_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .045, yend = .045 * systolicBP_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmHg)") +
  ylab("Standard Error (mmHg)") +
  ggtitle("Systolic BP") +
  xlim(0, .045) +
  ylim(0, 13)

testosterone_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "testosterone",
  n = 147477,
  m = 143937,
  n_samp = 100000
)

testosterone_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "testosterone",
  n = 147477,
  m = 143937,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_testosterone <- ggplot(data = testosterone_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .85, yend = .85 * testosterone_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (nmol/L)") +
  ylab("Standard Error (nmol/L)") +
  ggtitle("Testosterone") +
  xlim(0, .85) +
  ylim(0, 4)

wth_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "waist_to_hip",
  n = 180672,
  m = 155809,
  n_samp = 100000
)

wth_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "waist_to_hip",
  n = 180672,
  m = 155809,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_wth <- ggplot(data = wth_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .01, yend = .01 * wth_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects") +
  ylab("Standard Error") +
  ggtitle("Waist-to-Hip") +
  xlim(0, .01) +
  ylim(0, .06)

urate_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "urate",
  n = 148720,
  m = 172264,
  n_samp = 100000
)

urate_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "urate",
  n = 148720,
  m = 172264,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_urate <- ggplot(data = urate_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 12, yend = 12 * urate_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (umol/L)") +
  ylab("Standard Error (umol/L)") +
  ggtitle("Urate") +
  xlim(0, 12) +
  ylim(0, 60)

urea_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "urea",
  n = 148800,
  m = 172387,
  n_samp = 100000
)

urea_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "urea",
  n = 148800,
  m = 172387,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_urea <- ggplot(data = urea_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .0175, yend = .0175 * urea_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (mmol/L)") +
  ylab("Standard Error (mmol/L)") +
  ggtitle("Urea") +
  xlim(0, .0175) +
  ylim(0, 1.25)

waist_circ_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "waist_circ",
  n = 155850,
  m = 180701,
  n_samp = 100000
)

waist_circ_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "waist_circ",
  n = 155850,
  m = 180701,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_waist_circ <- ggplot(data = waist_circ_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .125, yend = .125 * waist_circ_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (cm)") +
  ylab("Standard Error (cm)") +
  ggtitle("Waist Circumference") +
  xlim(0, .125) +
  ylim(0, 10)

weight_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "weight",
  n = 155622,
  m = 180517,
  n_samp = 100000
)

weight_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "weight",
  n = 155622,
  m = 180517,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_weight <- ggplot(data = weight_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .06, yend = .06 * weight_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (Kg)") +
  ylab("Standard Error (Kg)") +
  ggtitle("Weight") +
  xlim(0, .06) +
  ylim(0, 9)

whole_body_fat_mass_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "whole_body_fat_mass",
  n = 152671,
  m = 178006,
  n_samp = 100000
)

whole_body_fat_mass_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "whole_body_fat_mass",
  n = 152671,
  m = 178006,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_whole_body_fat_mass <- ggplot(data = whole_body_fat_mass_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .05, yend = .05 * whole_body_fat_mass_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (Kg)") +
  ylab("Standard Error (Kg)") +
  ggtitle("Whole Body Fat Mass") +
  xlim(0, .05) +
  ylim(0, 7.5)

wth_bmi_adj_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "wth_bmi_adj",
  n = 155512,
  m = 180420,
  n_samp = 100000
)

wth_bmi_adj_info_top_hits <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "wth_bmi_adj",
  n = 155512,
  m = 180420,
  n_samp = 100000,
  pval_thresh = 5e-8
)

g_wth_bmi_adj <- ggplot(data = wth_bmi_adj_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .0025, yend = .0025 * wth_bmi_adj_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects") +
  ylab("Standard Error") +
  ggtitle("Waist-to-Hip BMI Adjusted") +
  xlim(0, .0025) +
  ylim(0, .1)

library(ggpubr)

fig <- ggarrange(
  g_albumin, g_arm_L, g_arm_R, g_bmi, 
  g_calcium, g_creatinine, g_diastolicBP, g_eosinophil, 
  g_fvc, g_hba1c, g_height, g_hip_circ,
  g_igf1, g_lymphocyte_perc, g_protein_total, g_pulse_rate, 
  g_rbc, g_shbg, g_systolicBP, g_testosterone, 
  g_urate, g_urea, g_waist_circ, g_weight,
  g_whole_body_fat_mass, g_wth, g_wth_bmi_adj,
  nrow = 7, ncol = 4
)

perc_better_vec <- 100 *c(
  albumin_info$pct_under_line,
  arm_fatfree_mass_L_info$pct_under_line,
  arm_fatfree_mass_R_info$pct_under_line,
  bmi_info$pct_under_line,
  calcium_info$pct_under_line,
  creatinine_info$pct_under_line,
  diastolicBP_info$pct_under_line,
  eosinophil_info$pct_under_line,
  fvc_info$pct_under_line,
  hba1c_info$pct_under_line,
  height_info$pct_under_line,
  hip_circ_info$pct_under_line,
  igf1_info$pct_under_line,
  lymphocyte_perc_info$pct_under_line,
  protein_total_info$pct_under_line,
  pulse_rate_info$pct_under_line,
  rbc_info$pct_under_line,
  shbg_info$pct_under_line,
  systolicBP_info$pct_under_line,
  testosterone_info$pct_under_line,
  urate_info$pct_under_line,
  urea_info$pct_under_line,
  waist_circ_info$pct_under_line,
  weight_info$pct_under_line,
  whole_body_fat_mass_info$pct_under_line,
  wth_info$pct_under_line,
  wth_bmi_adj_info$pct_under_line
)

names(perc_better_vec) <- c(
  "albumin",
  "arm fatfree mass left",
  "arm fatfree mass right",
  "bmi",
  "calcium",
  "creatinine",
  "diastolic BP",
  "eosinophil %",
  "fvc",
  "hba1c",
  "height",
  "hip circumference",
  "igf1",
  "lymphocyte %",
  "protein total",
  "pule rate",
  "red blood count",
  "shbg",
  "systolic BP",
  "testosterone",
  "urate",
  "urea",
  "waist circumference",
  "weight",
  "whole body fat mass",
  "waist to hip",
  "waist to hip bmi adjusted"
)

perc_better_vec <- sort(perc_better_vec)

perc_better_vec_top_hits <- 100 *c(
  albumin_info_top_hits$pct_under_line,
  arm_fatfree_mass_L_info_top_hits$pct_under_line,
  arm_fatfree_mass_R_info_top_hits$pct_under_line,
  bmi_info_top_hits$pct_under_line,
  calcium_info_top_hits$pct_under_line,
  creatinine_info_top_hits$pct_under_line,
  diastolicBP_info_top_hits$pct_under_line,
  eosinophil_info_top_hits$pct_under_line,
  fvc_info_top_hits$pct_under_line,
  hba1c_info_top_hits$pct_under_line,
  height_info_top_hits$pct_under_line,
  hip_circ_info_top_hits$pct_under_line,
  igf1_info_top_hits$pct_under_line,
  lymphocyte_perc_info_top_hits$pct_under_line,
  protein_total_info_top_hits$pct_under_line,
  pulse_rate_info_top_hits$pct_under_line,
  rbc_info_top_hits$pct_under_line,
  shbg_info_top_hits$pct_under_line,
  systolicBP_info_top_hits$pct_under_line,
  testosterone_info_top_hits$pct_under_line,
  urate_info_top_hits$pct_under_line,
  urea_info_top_hits$pct_under_line,
  waist_circ_info_top_hits$pct_under_line,
  weight_info_top_hits$pct_under_line,
  whole_body_fat_mass_info_top_hits$pct_under_line,
  wth_info_top_hits$pct_under_line,
  wth_bmi_adj_info_top_hits$pct_under_line
)

names(perc_better_vec_top_hits) <- c(
  "albumin",
  "arm fatfree mass left",
  "arm fatfree mass right",
  "bmi",
  "calcium",
  "creatinine",
  "diastolic BP",
  "eosinophil %",
  "fvc",
  "hba1c",
  "height",
  "hip circumference",
  "igf1",
  "lymphocyte %",
  "protein total",
  "pule rate",
  "red blood count",
  "shbg",
  "systolic BP",
  "testosterone",
  "urate",
  "urea",
  "waist circumference",
  "weight",
  "whole body fat mass",
  "waist to hip",
  "waist to hip bmi adjusted"
)

perc_better_vec_top_hits <- sort(perc_better_vec_top_hits)

perc_better_df <- data.frame(
  y = c(seq(1, length(perc_better_vec)), seq(1, length(perc_better_vec_top_hits))), 
  perc_better = c(perc_better_vec, perc_better_vec_top_hits), 
  criterion = c(rep("all", length(perc_better_vec)), rep("p < 5e-8", length(perc_better_vec_top_hits)))
)

ggplot(data = perc_better_df) +
  geom_point(aes(x = perc_better, y = y, color = criterion)) +
  scale_y_continuous(
    breaks = seq(1, length(perc_better_vec)),
    labels = names(perc_better_vec)
  ) +
  ylab("Trait") +
  xlab("Percent of Effects Better Estimated with GxE Model") +
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1.8,"cm")
  )

png(
  file = "~/Documents/paper_gwas_amplification/images/all_traits.png",
  width = 3500,
  height = 7000,
  res = 150
)
fig
dev.off()


# png(
#   file = "~/Documents/paper_gwas_amplification/images/waist_to_hip_cv.png",
#   width = 2000,
#   height = 1500,
#   res = 300
# )
g1 <- ggplot(data = test_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .15, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .003, yend = .003 * test_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects") +
  ylab("Standard Error") +
  ggtitle("Waist-to-Hip Ratio Decision Rule") +
  xlim(0, .003) +
  ylim(0, .04)
# dev.off()

test_info <- make_dec_rule_male_df(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "testosterone",
  n = 147477,
  m = 143937,
  n_samp = 100000
)

# now, get crossval plot

test_cv_info <- make_crossval_plot_male(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "testosterone"
)

g3 <- ggplot(data = test_cv_info, 
       aes(x = expected_advantage, y = unweighted_dec_advantage_add_male,
           ymin = unweighted_dec_advantage_add_male - sqrt(sd / samp_size),
           ymax = unweighted_dec_advantage_add_male + sqrt(sd / samp_size)
       )) +
  geom_point() +
  geom_errorbar() +
  xlab("Predicted MSE Advantage of Additive Model in Held-out Data") +
  ylab("MSE Advantage of Additive Model in Held-out Data") +
  ggtitle("Cross-Validation for Testosterone") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8))

library(ggplot2)
# 
# png(
#   file = "~/Documents/paper_gwas_amplification/images/testosterone_cv.png",
#   width = 2000,
#   height = 1500,
#   res = 300
# )
g2 <- ggplot(data = test_info$cv_df, aes(x = fx_diff, y = se)) +
  geom_point(alpha = .1, color = "orange", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = .2, yend = .2 * test_info$slope), linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8)) +
  xlab("Difference in Environment Specific Effects (ng/dl)") +
  ylab("Standard Error (ng/dl)") +
  ggtitle("Testosterone Decision Rule") +
  xlim(0, .2) +
  ylim(0, .75)
# dev.off()

wth_cv_info <- make_crossval_plot_female(
  data_dir = "/Users/ericweine/Documents/gwas_amplification/paper_figures/gxsex_cv_results",
  trait = "waist_to_hip"
)

g4 <- ggplot(data = wth_cv_info, 
       aes(x = expected_advantage, y = unweighted_dec_advantage_add_female,
           ymin = unweighted_dec_advantage_add_female - sqrt(sd / samp_size),
           ymax = unweighted_dec_advantage_add_female + sqrt(sd / samp_size)
       )) +
  geom_point() +
  geom_errorbar() +
  xlab("Predicted MSE Advantage of Additive Model in Held-out Data") +
  ylab("MSE Advantage of Additive Model in Held-out Data") +
  ggtitle("Cross-Validation for Waist-to-Hip Ratio") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 9), axis.title = element_text(size = 8))

library(ggpubr)
png(
  file = "~/Documents/paper_gwas_amplification/images/fig3.png",
  width = 5000,
  height = 5000,
  res = 300
)
fig <- ggarrange(
  g1, g2, g4, g3, nrow = 2, ncol = 2, labels = "AUTO"
)
annotate_figure(fig, top = textGrob("Application of Decision Rule to Sex-Stratified Human GWAS Data", hjust = .5, gp = gpar(cex = 1, fontsize = 25)))
dev.off()
