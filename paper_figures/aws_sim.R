set.seed(8675309)
'%>%' <- magrittr::'%>%'
library(doParallel)

Sigma <- list()
Sigma[[1]] <- matrix(data = 1, nrow = 2, ncol = 2)
Sigma[[2]] <- matrix(data = c(2.5, sqrt(2.5), sqrt(2.5), 1), ncol = 2)

maf_simulator <- function(n) {
  
  pmax(pmin(rbeta(n, .5, 10), .5 - 1e-05), 1e-5)
  
}

n.cores <- parallel::detectCores() 
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK",
  outfile = ""
)

doParallel::registerDoParallel(cl = my.cluster)

# metric_vec <- c()
# model_type_vec <- c()
# selection_method_vec <- c()
# value_vec <- c()
# amp_vec <- c()

amp_vec <- seq(0, 1, length.out = 8)

sim_out_list <- foreach::foreach(
  i = 1:length(amp_vec),
  .packages = c("gamp"),
  .verbose = T
) %do% {
  
  gamp::simulate_pgs_corr_fast_v2(
    n_e0 = 500000,
    n_e1 = 500000,
    n_snps = 10000,
    maf_simulator,
    e0_h2 = .25,
    e1_h2 = .25,
    Sigma = Sigma,
    pi = c(1 - amp_vec[i], amp_vec[i]),
    num_sims = 1,
    s = .5,
    pval_thresh_additive = .05,
    pval_thresh_GxE = .05,
    selection_methods = c("additive", "GxE")
  )
  
}

for (i in 1:length(amp_vec)) {
  
  readr::write_rds(sim_out_list[[i]], glue::glue("rds_data/polymetric_{amp_vec[i]}_zhu_hp_test.rds"))
  
}

# corr_vec <- c()
# amp_vec <- c()
# selection_method_vec <- c()
# beta_method_vec <- c()
# 
# for (selection_method in c("additive", "GxE")) {
#   
#   for (beta_method in c("additive", "GxE")) {
#     
#     for (amp in seq(0, 1, length.out = 8)) {
#       
#       sim_data <- readr::read_rds(
#         glue::glue("rds_data/polymetric_{amp}_zhu_hp_test.rds")
#       )
#       
#       amp_vec <- c(amp_vec, amp)
#       selection_method_vec <- c(selection_method_vec, selection_method)
#       beta_method_vec <- c(beta_method_vec, beta_method)
#       corr_vec <- c(corr_vec, sim_data[[selection_method]][[beta_method]][['corr']])
#       
#     }
#     
#   }
#   
# }
# 
# corr_df_hp <- data.frame(
#   selection_method = selection_method_vec,
#   beta_method = beta_method_vec,
#   amp = amp_vec,
#   corr = corr_vec
# )
# 
# corr_df_hp <- corr_df_hp %>%
#   dplyr::mutate(
#     selection_method = dplyr::if_else(
#       selection_method == "additive", 
#       "Additive Pval < .05",
#       "GxE Pval < .05"
#     ),
#     beta_method = dplyr::if_else(beta_method == "additive", "Additive", "GxE")
#   )
# 
# 
# hp_plot <- ggplot(data =corr_df_hp) +
#   geom_point(aes(x = amp, y = corr, color = selection_method), size = 1) +
#   geom_line(aes(x = amp, y = corr, color = selection_method, linetype = beta_method)) +
#   xlab("Percent of SNPs Amplified 2.5x in Environment 1") +
#   ylab("Test Set PGS Correlation") +
#   labs(color="Selection Criterion", linetype = "Estimate Model") +
#   ggtitle("PGS with h2 = .6 and n = 40,000 in Each Environment") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# corr_vec <- c()
# amp_vec <- c()
# selection_method_vec <- c()
# beta_method_vec <- c()
# 
# for (selection_method in c("additive", "GxE")) {
#   
#   for (beta_method in c("additive", "GxE")) {
#     
#     for (amp in seq(0, 1, length.out = 8)) {
#       
#       sim_data <- readr::read_rds(
#         glue::glue("rds_data/polymetric_{amp}_zhu_lp_test.rds")
#       )
#       
#       amp_vec <- c(amp_vec, amp)
#       selection_method_vec <- c(selection_method_vec, selection_method)
#       beta_method_vec <- c(beta_method_vec, beta_method)
#       corr_vec <- c(corr_vec, sim_data[[selection_method]][[beta_method]][['corr']])
#       
#     }
#     
#   }
#   
# }
# 
# corr_df_lp <- data.frame(
#   selection_method = selection_method_vec,
#   beta_method = beta_method_vec,
#   amp = amp_vec,
#   corr = corr_vec
# )
# 
# # Something here should get at the fact that there are two aspects to what's going on
# # one is selection and the other is estimation. 
# 
# corr_df_lp <- corr_df_lp %>%
#   dplyr::mutate(
#     selection_method = dplyr::if_else(
#       selection_method == "additive", 
#       "Additive Pval < .1",
#       "GxE Pval < .1"
#     ),
#     beta_method = dplyr::if_else(beta_method == "additive", "Additive", "GxE")
#   )
# 
# lp_plot <- ggplot(data =corr_df_lp) +
#   geom_point(aes(x = amp, y = corr, color = selection_method), size = 1) +
#   geom_line(aes(x = amp, y = corr, color = selection_method, linetype = beta_method)) +
#   xlab("Percent of SNPs Amplified 2.5x in Environment 1") +
#   ylab("Test Set PGS Correlation") +
#   labs(color="Selection Criterion", linetype = "Estimate Model") +
#   ggtitle("PGS with h2 = .6 and n = 1,000 in Each Environment") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggpubr::ggarrange(lp_plot, hp_plot)
