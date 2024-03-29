% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2ranges.R
\name{vcf2vranges}
\alias{vcf2vranges}
\title{Load a VCF file}
\usage{
vcf2vranges(
  vcfFile = NULL,
  genoField = NULL,
  translateBases = TRUE,
  genome = "hg38",
  phased = FALSE,
  region = NULL,
  sample = NULL
)
}
\arguments{
\item{vcfFile}{A path to a VCF file to be loaded.}

\item{genoField}{A vector of genotype IDs to be loaded from VCF [e.g. 'GT']}

\item{translateBases}{Set to \code{TRUE} if REF and ALT alleles should be reported as A,C,G or T.}

\item{genome}{A reference genome used by \link[VariantAnnotation]{readVcfAsVRanges} function. [e.g. 'hg38' - human]}

\item{phased}{If set to \code{TRUE} all unphased variants are removed.}

\item{region}{A \code{\link{GRanges}} object of genomic regions to be loaded from input VCF file.}

\item{sample}{A user defined set of sample IDs to be loaded.}
}
\value{
A \code{\link{VRanges-class}} object.
}
\description{
This function loads raw VCF file into a \code{\link{VRanges-class}} object.
}
\author{
David Porubsky
}
