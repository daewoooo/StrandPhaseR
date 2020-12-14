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
  
  ## Filter chromosomes to use	
  if (!is.null(chromosome)) { 	  
	if ( grepl(vcf$V1[1], pattern = "^chr|^cluster") & grepl(chromosome[1], pattern = "^chr|^cluster") ) {	
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

  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\/|\\|")
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

#' Load a VCF file
#' 
#' This function loads raw VCF file into a \code{\link{VRanges-class}} object.
#' 
#' @param vcfFile A path to a VCF file to be loaded.
#' @param genoField A vector of genotype IDs to be loaded from VCF [e.g. 'GT']
#' @param translateBases Set to \code{TRUE} if REF and ALT alleles should be reported as A,C,G or T.
#' @param genome A reference genome used by \link[VariantAnnotation]{readVcfAsVRanges} function. [e.g. 'hg38' - human]
#' @param phased If set to \code{TRUE} all unphased variants are removed.
#' @param region A \code{\link{GRanges}} object of genomic regions to be loaded from input VCF file.
#' @importFrom tidyr separate
#' @importFrom VariantAnnotation readVcfAsVRanges ScanVcfParam
#' @return A \code{\link{VRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
vcf2vranges <- function(vcfFile=NULL, genoField=NULL, translateBases=TRUE, genome='hg38', phased=FALSE, region=NULL) {
  ## Load vcf file into vranges object
  if (!is.null(region)) {
    if (all(is.character(genoField) & nchar(genoField) > 0)) {
      suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome, 
                                                                           param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA, geno = genoField, which = region)) )
    } else {
      suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome, 
                                                                           param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA, which = region)) )
    }
  } else {
    if (all(is.character(genoField) & nchar(genoField) > 0)) {
      suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome, 
                                                                           param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA, geno = genoField)) )
    } else {
      suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome, 
                                                                           param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA)) )
    }
  }  
  ## Remove unhased variants
  if (phased) {
    mask <- grepl(vcf.vranges$GT, pattern = "\\/")
    vcf.vranges <- vcf.vranges[!mask]
  }
  ## Split genotype field into hapltypes
  gen.field <-  tidyr::separate(data = as(mcols(vcf.vranges), 'data.frame'), col = GT, into = c("H1", "H2"))
  ## Translate 0/1 haplotypes into nucleotides
  if (translateBases) {
    allele1 <- rep(".", length(vcf.vranges))
    allele2 <- rep(".", length(vcf.vranges))
    allele1[gen.field$H1 == 0] <- ref(vcf.vranges)[gen.field$H1 == 0]
    allele1[gen.field$H1 == 1] <- alt(vcf.vranges)[gen.field$H1 == 1]
    allele1[allele1 == "."] <- 'N'
    allele2[gen.field$H2 == 0] <- ref(vcf.vranges)[gen.field$H2 == 0]
    allele2[gen.field$H2 == 1] <- alt(vcf.vranges)[gen.field$H2 == 1]
    allele2[allele2 == "."] <- 'N'
    gen.field$H1 <- allele1
    gen.field$H2 <- allele2
  }
  
  if (!phased) {
    colnames(gen.field) <- c('allele1', 'allele2')
  }
  
  ## Construct final VRanges object
  mcols(vcf.vranges) <- gen.field
  ## Keep only seqlevels present in the final VRanges object
  vcf.vranges <- keepSeqlevels(vcf.vranges, value = runValue(seqnames(vcf.vranges)))
  
  return(vcf.vranges)
} 
