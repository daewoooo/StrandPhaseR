#' Read VCF file into a \code{\link{GRanges}} object
#' 
#' @param vcfFile VCF file containing phased alleles
#' @param genotypeField In case of multiple samples phased in a single file (Default: 1)
#' @inheritParams phaseChromosome 
#' 
#' @author David Porubsky
#' @export

vcf2ranges <- function(vcfFile=NULL, genotypeField=1, chromosome=NULL) {
  message("Loading VCF file ...", appendLF=F); ptm <- proc.time()
  
  vcf <- read.table(vcfFile, stringsAsFactors = FALSE, fill=TRUE)
  
  #filter chromosomes to use	
  if (!is.null(chromosome)) { 	  
	if ( grepl(vcf$V1[1], pattern = "^chr") & grepl(chromosome[1], pattern = "^chr") ) {	
		vcf <- vcf[vcf$V1 %in% chromosome,]
	} else {
		stop('Specified chromosomes names do not match VCF chromosome names!!!')
  	}
  }	

  #remove indels
  ref.len <- sapply(vcf$V4, nchar, USE.NAMES = F)
  alt.len <- sapply(vcf$V5, nchar, USE.NAMES = F)
  vcf <- vcf[ref.len == 1 & alt.len == 1,]
  
  column <- genotypeField + 9
  vcf <- vcf[,c(1:9, column)]

  #gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- lapply(as.character(vcf[,10]), function(x) strsplit(x,':')[[1]][1])
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block,"\\/|\\|")
  alleles <- do.call(rbind, alleles)

  #filter only HET positions
  mask <- alleles[,1] != alleles[,2]
  #gen.block <- gen.block[mask,]
  vcf <- vcf[mask,]
  
  #export vcf in GRanges object
  vcf.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=as.numeric(vcf$V2), end=as.numeric(vcf$V2)), ref=vcf$V4, alt=vcf$V5)
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  return(vcf.gr)
}