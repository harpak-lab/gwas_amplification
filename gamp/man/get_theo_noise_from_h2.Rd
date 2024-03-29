% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygenic.R
\name{get_theo_noise_from_h2}
\alias{get_theo_noise_from_h2}
\title{Get Observation Noise Based on Heritability}
\usage{
get_theo_noise_from_h2(maf_vec, fx_vec, h2)
}
\arguments{
\item{maf_vec}{vector containing the minor allele frequency of each snp}

\item{fx_vec}{vector containing true effect size in each snp}

\item{h2}{number indiciating desired heritability}
}
\value{
environmental variance
}
\description{
Get Observation Noise Based on Heritability
}
