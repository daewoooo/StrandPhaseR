#' Read VCF file into a \code{\link{GRanges}} object
#' 
#' @param vcfFile VCF file containing phased alleles
#' @param genotypeField In case of multiple samples phased in a single file (Default: 1)
#' 
#' @author David Porubsky
#' @export

vcf2ranges <- function(vcfFile=NULL, genotypeField=1, chromosomes=NULL) {
  vcf <- read.table(vcfFile, stringsAsFactors = F, fill=TRUE)
  
  #filter chromosomes to use	
  if (!is.null(chromosomes)) { 	  
	if ( grepl(vcf$V1[1], pattern = "^chr") & grepl(chromosomes[1], pattern = "^chr") ) {	
		vcf <- vcf[vcf$V1 %in% chromosomes,]
	} else {
		stop('Specified chromosomes names do not match VCF chromosome names!!!')
  	}
  }	

  column <- genotypeField + 9
  vcf <- vcf[,c(1:9, column)]
  
  #remove indels
  ref.len <- sapply(vcf$V4, nchar, USE.NAMES = F)
  alt.len <- sapply(vcf$V5, nchar, USE.NAMES = F)
  vcf <- vcf[ref.len == 1 & alt.len == 1,]	

  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\/")
  alleles <- do.call(rbind, alleles)

  #filter only HET positions
  mask <- alleles[,1] != alleles[,2]  
  vcf <- vcf[mask, c(1,2)]
  
  vcf.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=as.numeric(vcf$V2), end=as.numeric(vcf$V2)))
  
  return(vcf.gr)
}
