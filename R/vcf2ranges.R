#' Read VCF file into a \code{\link{GRanges}} object
#' 
#' @param vcfFile VCF file containing phased alleles
#' @param genotypeField In case of multiple samples phased in a single file (Default: 1)
#' 
#' @author David Porubsky
#' @export

vcf2ranges <- function(vcfFile=NULL, genotypeField=1) {
  vcf <- read.table(vcfFile, stringsAsFactors = F)
  
  column <- genotypeField + 9
  vcf <- vcf[,c(1:9, column)]
  
  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\/")
  alleles <- do.call(rbind, alleles)

  #filter only HET positions
  mask <- alleles[,1] != alleles[,2]  
  vcf <- vcf[mask, c(1,2)]
  
  vcf.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2))
  
  return(vcf.gr)
}
