% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bamregion2GRanges.R
\name{bamregion2GRanges}
\alias{bamregion2GRanges}
\title{Import BAM file into GRanges}
\usage{
bamregion2GRanges(
  bamfile,
  bamindex = bamfile,
  region = NULL,
  pairedEndReads = FALSE,
  min.mapq = 10,
  filterAltAlign = TRUE
)
}
\arguments{
\item{bamindex}{Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.}

\item{region}{If only a subset of the genomic regions should be loaded.}

\item{pairedEndReads}{Set to \code{TRUE} if you have paired-end reads in your file.}

\item{min.mapq}{Minimum mapping quality when importing from BAM files.}

\item{file}{Bamfile with aligned reads.}
}
\description{
Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
}
\author{
Aaron Taudt, David Porubsky, Ashley Sanders
}
