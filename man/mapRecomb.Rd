% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapRecomb.R
\name{mapRecomb}
\alias{mapRecomb}
\title{Map meiotic recombination events in a single family trio [father-mother-child].}
\usage{
mapRecomb(
  parent1 = NULL,
  parent2 = NULL,
  child = NULL,
  genome = "hg38",
  method = "CBS",
  minSeg = 100,
  smooth = 3,
  collapse.amb = TRUE
)
}
\arguments{
\item{parent1}{A path to a VCF file to be loaded for a parent 1.}

\item{parent2}{A path to a VCF file to be loaded for a parent 2.}

\item{child}{A path to a VCF file to be loaded for a child.}

\item{genome}{A reference genome used by \link[VariantAnnotation]{readVcfAsVRanges} function. [e.g. 'hg38' - human]}

\item{method}{A user defined method to be used to map changes in haplotype blocks [default: CBS]}

\item{minSeg}{Minimal length (number of variants) being reported as haplotype block (\code{fastseg} parameter).}

\item{smooth}{Number of consecutive variants being considered as a random error and so being corrected (flipped).}

\item{collapse.amb}{Set to \code{TRUE} if segments with ambiguous haplotype assignments should be collapsed.}
}
\value{
A \code{list} object that contains mapped meiotic breakpoints and inherited hapltype segments 
for each homolog in a child.
}
\description{
This function takes as an input phased VCF files for each member of a family trio and maps breakpoints
of meiotic recombination in each inherited parental homolog.
}
\author{
David Porubsky
}
