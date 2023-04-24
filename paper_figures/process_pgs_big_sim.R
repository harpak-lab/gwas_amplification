# here, I want to process the polygenic score simulations
# first, I want to recall what one object looks like

num_replicates <- 22
amp_perc_vec <- c()
replicate_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()
corr_vec <- c()

for (amp_perc in seq(0, 1, .1)) {
  
  for (replicate in 1:num_replicates) {
    
    sim <- readr::read_rds(
      glue::glue("pgs_sim_res/pgs_sim_n_1000_n_snps_5000_h2_0.4_amp_coef_1.5",
                 "_amp_perc_{amp_perc}_n_sims_1300_replicate_{replicate}_fixed_asc_snps.rds")
    )
    
    for (selection_method in c("additive", "GxE")) {
      
      for (beta_method in c("additive", "GxE")) {
        
        amp_perc_vec <- c(amp_perc_vec, amp_perc)
        replicate_vec <- c(replicate_vec, replicate)
        selection_method_vec <- c(selection_method_vec, selection_method)
        beta_method_vec <- c(beta_method_vec, beta_method)
        corr_vec <- c(corr_vec, sim[[selection_method]][[beta_method]][["corr"]])
        
      }
      
    }
    
  }
  
}

sim_df_low_p <- data.frame(
  amp_perc = amp_perc_vec,
  replicate = replicate_vec,
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  corr = corr_vec
)
  
# here, I want to process the polygenic score simulations
# first, I want to recall what one object looks like

num_replicates <- 22
amp_perc_vec <- c()
replicate_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()
corr_vec <- c()

for (amp_perc in seq(0, 1, .1)) {
  
  for (replicate in 1:num_replicates) {
    
    sim <- readr::read_rds(
      glue::glue("pgs_sim_res/pgs_sim_n_50000_n_snps_5000_h2_0.4_amp_coef_1.5",
                 "_amp_perc_{amp_perc}_n_sims_1300_replicate_{replicate}_fixed_asc_snps.rds")
    )
    
    for (selection_method in c("additive", "GxE")) {
      
      for (beta_method in c("additive", "GxE")) {
        
        amp_perc_vec <- c(amp_perc_vec, amp_perc)
        replicate_vec <- c(replicate_vec, replicate)
        selection_method_vec <- c(selection_method_vec, selection_method)
        beta_method_vec <- c(beta_method_vec, beta_method)
        corr_vec <- c(corr_vec, sim[[selection_method]][[beta_method]][["corr"]])
        
      }
      
    }
    
  }
  
}

sim_df_high_p <- data.frame(
  amp_perc = amp_perc_vec,
  replicate = replicate_vec,
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  corr = corr_vec
)

# now, get the mash results

num_replicates <- 22
amp_perc_vec <- c()
replicate_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()
corr_vec <- c()

for (amp_perc in seq(0, 1, .1)) {
  
  for (replicate in 1:num_replicates) {
      
    if (!(replicate == 16 && amp_perc == 0.4) && !(replicate == 24 && amp_perc == 0) &&
        !(replicate == 26 && amp_perc == 0) && !(replicate == 29 && amp_perc == 0) &&
        !(replicate == 31 && amp_perc == 0) && !(replicate == 32 && amp_perc == 0) &&
        !(replicate == 33 && amp_perc == 0) && !(replicate == 24 && amp_perc == 0.1) &&
        !(replicate == 30 && amp_perc == 0.1) && !(replicate == 31 && amp_perc == 0.1) &&
        !(replicate == 25 && amp_perc == 0.2) && !(replicate == 26 && amp_perc == 0.2)
        ){
    
    sim <- readr::read_rds(
        glue::glue("pgs_sim_res/pgs_sim_n_1000_n_snps_5000_h2_0.4_amp_coef_1.5",
                   "_amp_perc_{amp_perc}_n_sims_100_replicate_{replicate}_fixed_asc_snps_mash.rds")
      )
      
    
    for (selection_method in c("mash")) {
      
      for (beta_method in c("mash")) {
        
        amp_perc_vec <- c(amp_perc_vec, amp_perc)
        replicate_vec <- c(replicate_vec, replicate)
        selection_method_vec <- c(selection_method_vec, selection_method)
        beta_method_vec <- c(beta_method_vec, beta_method)
        corr_vec <- c(corr_vec, sim[[selection_method]][[beta_method]][["corr"]])
        
      }
      
    }
    
    }
    
  }
  
}

sim_df_low_p_mash <- data.frame(
  amp_perc = amp_perc_vec,
  replicate = replicate_vec,
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  corr = corr_vec
)

num_replicates <- 22
amp_perc_vec <- c()
replicate_vec <- c()
selection_method_vec <- c()
beta_method_vec <- c()
corr_vec <- c()

for (amp_perc in seq(0, 1, .1)) {
  
  for (replicate in 1:num_replicates) {
    
    if (!(replicate == 7 && amp_perc == 0.5)  && !(replicate == 3 && amp_perc == 1)) {
    
    sim <- readr::read_rds(
      glue::glue("pgs_sim_res/pgs_sim_n_50000_n_snps_5000_h2_0.4_amp_coef_1.5",
                 "_amp_perc_{amp_perc}_n_sims_100_replicate_{replicate}_fixed_asc_snps_mash.rds")
    )
    
    for (selection_method in c("mash")) {
      
      for (beta_method in c("mash")) {
        
        amp_perc_vec <- c(amp_perc_vec, amp_perc)
        replicate_vec <- c(replicate_vec, replicate)
        selection_method_vec <- c(selection_method_vec, selection_method)
        beta_method_vec <- c(beta_method_vec, beta_method)
        corr_vec <- c(corr_vec, sim[[selection_method]][[beta_method]][["corr"]])
        
      }
      
    }
    
    }
    
  }
  
}

sim_df_high_p_mash <- data.frame(
  amp_perc = amp_perc_vec,
  replicate = replicate_vec,
  selection_method = selection_method_vec,
  beta_method = beta_method_vec,
  corr = corr_vec
)

# the correct files are not loaded in above 
# but once they are I can make the plot
# which will hopefully happen soon
library(dplyr)
corr_df_low_p <- sim_df_low_p %>%
  rbind(sim_df_low_p_mash) %>%
  dplyr::group_by(amp_perc, selection_method, beta_method) %>%
  dplyr::summarize(
    mean_corr = mean(corr)
  )

corr_df_high_p <- sim_df_high_p %>%
  rbind(sim_df_high_p_mash) %>%
  dplyr::group_by(amp_perc, selection_method, beta_method) %>%
  dplyr::summarize(
    mean_corr = mean(corr)
  )

corr_df_high_p <- corr_df_high_p %>%
  dplyr::mutate(
    selection_method = dplyr::case_when(
      selection_method == "additive" ~ "additive p-value",
      selection_method == "GxE" ~ "GxE p-value",
      selection_method == "mash" ~ "polygenic adaptive shrinkage\n pseudo p-value"
    )) %>%
  dplyr::mutate(
    beta_method = if_else(
      beta_method == "mash", 
      "polygenic adapditve shrinkage",
      beta_method
    )
  ) %>%
  dplyr::filter(!(beta_method == "additive" & selection_method == "GxE p-value")) %>%
  dplyr::mutate(comb_method = paste(beta_method, selection_method))

corr_df_low_p <- corr_df_low_p %>%
  dplyr::mutate(
    selection_method = dplyr::case_when(
      selection_method == "additive" ~ "additive p-value",
      selection_method == "GxE" ~ "GxE p-value",
      selection_method == "mash" ~ "polygenic adaptive shrinkage\npseudo p-value"
    )) %>%
  dplyr::mutate(
    beta_method = if_else(
      beta_method == "mash", 
      "polygenic adaptive shrinkage",
      beta_method
    )
  ) %>%
  dplyr::filter(!(beta_method == "additive" & selection_method == "GxE p-value")) %>%
  dplyr::mutate(comb_method = paste(beta_method, selection_method))

library(ggplot2)

g1 <- ggplot(data = corr_df_low_p) +
  geom_line(aes(x = amp_perc, y = mean_corr, color = comb_method), size = 1.25) +
  xlab("Percent of SNPs Amplified 1.5x") +
  ylab("Test Set PGS Correlation") +
  ggtitle("Low Power") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("blue", "magenta", "red", "green3")) +
  annotate("text", x = .75, y = .275, label = "polygenic adaptive shrinkage", color = "green3", size = 3.5) +
  geom_segment(aes(x = .75, y = .28, xend = .75, yend = .2995),
               arrow = arrow(length = unit(0.5, "cm")), color = "green3") +
  annotate("text", x = .85, y = .18, label = "GxE", color = "red", size = 3.5) +
  geom_segment(aes(x = .85, y = .175, xend = .85, yend = .1575),
               arrow = arrow(length = unit(0.5, "cm")), color = "red") +
  annotate("text", x = .35, y = .165, label = "GxE; ascertainment based on additive model", color = "magenta", size = 3.5) +
  geom_segment(aes(x = .35, y = .168, xend = .35, yend = .1925),
               arrow = arrow(length = unit(0.5, "cm")), color = "magenta") +
  annotate("text", x = .175, y = .255, label = "additive", color = "blue", size = 3.5) +
  geom_segment(aes(x = .175, y = .25, xend = .175, yend = .225),
               arrow = arrow(length = unit(0.5, "cm")), color = "blue")

g2 <- ggplot(data = corr_df_high_p) +
  geom_line(aes(x = amp_perc, y = mean_corr, color = comb_method), size = 1.25) +
  xlab("Percent of SNPs Amplified 1.5x") +
  ylab("Test Set PGS Correlation") +
  ggtitle("High Power") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("blue", "magenta", "red", "green3")) +
  annotate("text", x = .8, y = .595, label = "polygenic adaptive shrinkage", color = "green3", size = 3.5) +
  geom_segment(aes(x = .8, y = .5955, xend = .8, yend = .60275),
               arrow = arrow(length = unit(0.5, "cm")), color = "green3") +
  annotate("text", x = .075, y = .592, label = "additive", color = "blue", size = 3.5) +
  geom_segment(aes(x = .075, y = .593, xend = .075, yend = .599),
               arrow = arrow(length = unit(0.5, "cm")), color = "blue") +
  annotate("text", x = .05, y = .58, label = "GxE", color = "red", size = 3.5) +
  geom_segment(aes(x = .05, y = .579, xend = .05, yend = .572),
               arrow = arrow(length = unit(0.5, "cm")), color = "red") +
  annotate("text", x = .53, y = .581, label = "GxE; ascertainment based on additive model", color = "magenta", size = 3.5) +
  geom_segment(aes(x = .53, y = .582, xend = .53, yend = .5905),
               arrow = arrow(length = unit(0.5, "cm")), color = "magenta")
  

library(ggpubr)
library(grid)

fig4 <- ggarrange(
  g1 + rremove("ylab") + rremove("xlab"),
  g2 + rremove("ylab") + rremove("xlab"), 
  nrow = 1,
  legend = "none",
  labels = "AUTO"
)

annotate_figure(fig4, 
                left = textGrob("Test Set PGS Correlation", rot = 90, vjust = 1, gp = gpar(cex = 1, fontsize = 12)),
                bottom = textGrob("Percent of Variant Effects Amplified 1.5x", gp = gpar(cex = 1, fontsize = 12), hjust = .45))



  