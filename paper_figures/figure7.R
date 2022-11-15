# here, want to examine distribution of effect sizes by MAF
# I don't think this should be too difficult
set.seed(1)
'%>%' <- magrittr::'%>%'
source("./utils.R")

male_df <- readr::read_tsv(
  file = '../data/male_all.height.glm.linear',
  col_select = c("ID", "BETA", "SE", "OBS_CT", "P", "#CHROM", "POS")
) %>%
  dplyr::rename(
    b_m = BETA, se_m = SE, n_m = OBS_CT, p_m = P, chrom = "#CHROM", pos = POS
  ) %>%
  dplyr::filter(p_m < 5e-8)

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
    dplyr::filter(ID %in% male_df$ID)

  all_c_df <- rbind(all_c_df, c_df)

}

male_df <- male_df %>%
  dplyr::inner_join(all_c_df, by = "ID")


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

out <- sample_snps(male_df, ld_df, samps_per_block = 25)

maf_low <- out %>%
  dplyr::filter(maf >= 0 & maf <= .5)

plot(density(maf_low$b_m), main = "Density of Sig Effects")

