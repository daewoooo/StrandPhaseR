% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correctInversionPhasing.R
\name{correctHetInv}
\alias{correctHetInv}
\title{Correct heterozygous inversion in StrandPhaseR VCF file.}
\usage{
correctHetInv(correct.gr = NULL, vcf.file = NULL, het.haps = NULL, ID = "")
}
\arguments{
\item{correct.gr}{A \code{\link{GRanges}} object of homozygous inversion site to be corrected.}

\item{vcf.file}{A \pkg{\link{StrandPhaseR}} formatted VCF file to be corrected for inversion phasing.}

\item{het.haps}{A \code{data.frame} with phased alleles per inverted and reference haplotype.}

\item{ID}{A unique id to be appended at the end of each corrected VCF file.}
}
\value{
\code{NULL}
}
\description{
This function takes as an input region deemed to be a heterozygous inversion along with phased SNVs within a region and uses
this information to add/correct into the \pkg{\link{StrandPhaseR}} formatted VCF file.
}
\author{
David Porubsky
}
