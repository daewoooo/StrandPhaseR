% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembleHaps.R
\name{assembleHaps}
\alias{assembleHaps}
\title{This will take sorted matrices and will calculate concordance of each single cell to the consensus haplotypes in order to assemble highly 
accurate haplotypes}
\usage{
assembleHaps(data.object, translateBases = FALSE, concordance = 0.9)
}
\arguments{
\item{data.object}{containing sorted watson and crick haplotypes of each single cell}

\item{translateBases}{translates integer coded bases (1,2,3,4) into letters (A,C,G,T)}

\item{concordance}{Level of agreement between single cell and consensus haplotypes}
}
\description{
This will take sorted matrices and will calculate concordance of each single cell to the consensus haplotypes in order to assemble highly 
accurate haplotypes
}
\author{
David Porubsky
}
