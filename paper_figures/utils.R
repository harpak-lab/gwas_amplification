'%>%' <- magrittr::'%>%'
library(data.table)

sample_snps <- function(snp_df, ld_df, samps_per_block = 1) {

  snp_df <- setDT(snp_df)
  ld_df <- setDT(ld_df)

  out <- snp_df[ld_df, on = .(chrom == chrom, pos >= start, pos <= stop),
                 nomatch = 0, ]

  samp <- out %>%
    dplyr::group_by(ld_block) %>%
    dplyr::sample_n(samps_per_block, replace = T) %>%
    dplyr::select(
      -pos, -pos.1, -ld_block
    ) %>%
    dplyr::distinct()

  return(samp)

}

