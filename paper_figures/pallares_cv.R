library(data.table)
library(aod)
library(glue)
library(dplyr)

#command_args = commandArgs(trailingOnly = TRUE)
#start = as.integer(command_args[1])
#stop = as.integer(command_args[2])
start <- 11501
stop <- 11750
print(start)
print(stop)

set.seed(1)
#data_file_prefix <- "/scratch/09251/eweine/pallares_cv/data"
data_file_prefix <- "/Users/ericweine/Downloads/doi_10.5061_dryad.hhmgqnkgr__v7"

# count data for beta binomial model
alt <- fread(glue('{data_file_prefix}/29Jun20_merged_alt_counts_allCHR.txt'))
both <- fread(glue('{data_file_prefix}/29Jun20_merged_both_counts_allCHR.txt'))
info <- read.delim(glue('{data_file_prefix}/29Jun20_merged_sample_info.txt'))

summary_table <- read.delim(glue('{data_file_prefix}/SummaryTable_allsites_12Nov20.txt'))

sites_df <- data.frame(stringr::str_split_fixed(summary_table$site, ":", 2))
rm("summary_table")
gc(verbose = FALSE)

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
)[start:stop]

rm("sites_df")
rm("sites_sample_df")
gc(verbose = FALSE)

# Now, I can do cross-validation just on the sites that are randomly selected
# I have 12K to work with, which seems reasonable
alt <- subset(alt, site %in% selected_sites)
both <- subset(both, site %in% selected_sites)

K <- 5
fold_assignments <- sample(1:K, nrow(info), replace = TRUE)
info$fold <- fold_assignments

info <- info %>%
  dplyr::mutate(
    time = case_when(timepoint == "TN" ~ 1,
                     timepoint == "T0" ~ 0)
    )

out_mat <- matrix(nrow = nrow(alt) * K, ncol = 10)
colnames(out_mat) <- c(
  "site_id", "fold", "b_t_train", "b_c_train", "b_hs_train", "se_t_train", 
  "se_c_train", "se_hs_train", "b_c_test", "b_hs_test"
)

# fit model and save pvals, coefs, se, lik at each of the sites
for (i in 1:nrow(alt)) {
  
  #browser()
  info$tmp_x <- as.vector(t(alt[i,-1]))
  info$tmp_y <- as.vector(t(both[i,-1]))
  
  out_vec_list <- list()
  print(i)
  
  if (i == 160) {
    
    print(0)
    
  }
  
  for (k in 1:5) {
    
    #browser()
    # exclude samples with no reads at a given site, and samples that couldn't be sexed
    tmp_info_train <- subset(info, tmp_y>0 & sex!='unknown' & fold != k)
    tmp_info_test <- subset(info, tmp_y>0 & sex!='unknown' & fold == k)
    
    additive_mod_train <- tryCatch(
      betabin(
        cbind(tmp_x, tmp_y - tmp_x) ~ time + sequencing_batch + meta_cage + sex,
        ~1, 
        data=tmp_info_train,
        control = list(maxit = 2500, reltol = 1e-7)
      ),
      error=function(x){}
    )
    
    GxE_mod_train <- tryCatch(
      betabin(
        cbind(tmp_x, tmp_y - tmp_x) ~ condition + sequencing_batch + meta_cage + sex,
        ~1, 
        data=tmp_info_train,
        control = list(maxit = 2500, reltol = 1e-7)
        ),
      error=function(x){}
    )
    
    GxE_mod_test <- tryCatch(
      betabin(
        cbind(tmp_x, tmp_y - tmp_x) ~ condition + sequencing_batch + meta_cage + sex,
        ~1, 
        data=tmp_info_test,
        control = list(maxit = 2500, reltol = 1e-7)
      ),
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"site_id"] <- i
    out_mat[K * (i - 1) + k,"fold"] <- k
    out_mat[K * (i - 1) + k,"b_t_train"] <- tryCatch(
      attributes(summary(additive_mod_train))$Coef["time", "Estimate"], 
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"b_hs_train"] <- tryCatch(
      attributes(summary(GxE_mod_train))$Coef["conditionHS", "Estimate"], 
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"b_c_train"] <- tryCatch(
      attributes(summary(GxE_mod_train))$Coef["conditionC", "Estimate"], 
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"se_t_train"] <- tryCatch(
      attributes(summary(additive_mod_train))$Coef["time", "Std. Error"],
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"se_hs_train"] <- tryCatch(
      attributes(summary(GxE_mod_train))$Coef["conditionHS", "Std. Error"], 
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"se_c_train"] <- tryCatch(
      attributes(summary(GxE_mod_train))$Coef["conditionC", "Std. Error"],
      error=function(x){}
    )
    
    out_mat[K * (i - 1) + k,"b_c_test"] <- tryCatch(
      attributes(summary(GxE_mod_test))$Coef["conditionC", "Estimate"], 
      error=function(x){
        print(glue("error on interation {i} and fold {k}"))
      }
    )
    
    out_mat[K * (i - 1) + k,"b_hs_test"] <- tryCatch(
      attributes(summary(GxE_mod_test))$Coef["conditionHS", "Estimate"], 
      error=function(x){}
    )
    
  }

}

out_df <- as.data.frame(out_mat)
out_df <- out_df %>%
  mutate(site = selected_sites[site_id]) %>%
  select(-site_id)

#readr::write_csv(out_mat, "/scratch/09251/eweine/pallares_cv/results/cv_df.csv")
readr::write_csv(
  out_df, 
  glue("~/Documents/gwas_amplification/paper_figures/pallares_cv_start_{start}_stop_{stop}.csv")
)
