#' Genotype Strand-seq cell based on population genotypes.
#'
#' This function will import variants (SNVs) covered in individual Strand-seq libraries in a form of BAM files aligned to a
#' a reference genome such as GRCh38. Such SNVs are compared to a population genotypes such as 1000G sample panel.
#' 
#' @param inputfolder A data folder where individual Strand-seq BAM files are stored.
#' @param strandS.vcf A VCF file that contains SNVs called from merged Strand-seq libraries using RTG tool.
#' @param popul.vcf.list A names list of paths to VCF files (per chromosome) that contains SNVs from multiple individuals such as 1000G sample panel.
#' @param wc.regions A Watson-crick regions per library and per chromosome as defined by breakpointR.
#' @param chromosomes List of chromosomes to be used for genotyping.
#' @param min.snv.cov A minumum number of Strand-seq reads required to cover a SNV position in 'strandS.vcf'.  
#' @param max.snv.cov A maximum number of Strand-seq reads allowed to cover a SNV position in 'strandS.vcf'.  
#' @param max.snv.per.chr A maximum number of SNVs to be loaded from population panel VCF ('popul.vcf) per chromosome.
#' @param blacklist A\code{\link{GRanges-class}} object of genomic regions to filter out SNV positions.
#' @return A \code{data.frame} object.
#' @importFrom S4Vectors endoapply
#' @importFrom dplyr filter group_by summarise n
#' @importFrom VariantAnnotation readVcfAsVRanges sampleNames ScanVcfParam ref alt
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom Rsamtools indexTabix
#' @importFrom stringr str_split_fixed
#' @author David Porubsky
#' @export
#' 
genotypeStrandScells <- function(inputfolder=NULL, strandS.vcf=NULL, popul.vcf.list=NULL, wc.regions=NULL, chromosomes=paste0('chr', c(1:22)), min.snv.cov=5, max.snv.cov=30, max.snv.per.chr=30000, blacklist=NULL) {
  ## Check user input
  if (any(is.null(inputfolder), is.null(strandS.vcf), is.null(popul.vcf.list), is.null(wc.regions))) {
    warning("Required parameters: 'inputfolder, 'strandS.vcf', 'popul.vcf.list', 'wc.regions' are not fully defined, aborting!!!")
  }
  
  ## Load and filter SNVs detected by RTG tool in merged Strand-seq cells
  message("Loading SNVs found in merged Strand-seq data ...")
  snvs <- StrandPhaseR::vcf2vranges(vcfFile = strandS.vcf, genoField = c('GT', 'DP'), translateBases = FALSE, phased = FALSE)
  ## Filter non-alt alleles
  #snvs <- snvs[!is.na(alt(snvs))]
  ## Keep HETs only
  snvs <- snvs[which(snvs$H1 != snvs$H2)]
  ## Filter sites based depth of coverage
  #snvs <- snvs[totalDepth(snvs) >= 5]
  DP <- VariantAnnotation::totalDepth(snvs)
  ## Remove seqlength and genome ID
  seqlengths(snvs) <- NA
  genome(snvs) <- NA
  ## Remove mutlinucleotide alts
  snvs <- snvs[nchar(VariantAnnotation::alt(snvs)) == 1]
  ## Convert to GRanges object
  snvs <- as(snvs, "GRanges")
  snvs$DP <- DP
  ## Remove any remaining mutlinucleotide positions
  snvs <- snvs[width(snvs) == 1]
  ## Remove variants from regions defined in blacklist
  if (!is.null(blacklist) & class(blacklist) == "GRanges") {
    snvs <- suppressWarnings( IRanges::subsetByOverlaps(snvs, blacklist, invert = TRUE) )
  }  
  
  #chromosomes <- 'chr1'
  compare.per.chr <- list()
  for (chr in chromosomes) {
    message("Working on chromosome: ", chr)
    ## Select SNVs for a single chromosome
    snvs.chr <- GenomeInfoDb::keepSeqlevels(snvs, value = chr, pruning.mode = 'coarse')
    ## Filter SNVs by depth
    if (all(!is.na(snvs.chr$DP))) {
      if (min.snv.cov > 0) {
        snvs.chr <- snvs.chr[snvs.chr$DP >= min.snv.cov]
      } 
      if (max.snv.cov > 0) {
        snvs.chr <- snvs.chr[snvs.chr$DP <= max.snv.cov]
      }
    }
    ## Filter SNVs by max allowed SNVs per chromosome
    if (max.snv.per.chr > 0) {
      if (length(snvs.chr) > max.snv.per.chr) {
        sub.pos <- sample(1:length(snvs.chr))[1:max.snv.per.chr]
        snvs.chr <- sort(snvs.chr[sub.pos])
      }  
    }
    
    ## Load WC regions
    WC.regions <- utils::read.table(wc.regions, header=FALSE)
    WC.regions <- GenomicRanges::GRanges(seqnames=WC.regions$V1, IRanges(start=WC.regions$V2, end=WC.regions$V3), filename=as.character(WC.regions$V4))
    WC.regions.chr <- GenomeInfoDb::keepSeqlevels(WC.regions, value = chr, pruning.mode = 'coarse')
    ## Keep one WC region per chromosome and per cell
    WC.regions.chr.l <- split(WC.regions.chr, WC.regions.chr$filename)
    WC.regions.chr.l <- S4Vectors::endoapply(WC.regions.chr.l, function(x) x[which.max(width(x))])
    WC.regions.chr <- unlist(WC.regions.chr.l)
    
    ## Extract haplotypes per WC regions
    matrices <- StrandPhaseR::loadMatrices(inputfolder=inputfolder, positions=snvs.chr, WCregions=WC.regions.chr, pairedEndReads=TRUE, min.mapq=10, min.baseq=20)
    
    message("    Loading 1000G SNVs ...")
    ## Select genomic sites to load SNVs
    target.sites <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges::IRanges(start=matrices$genomic.pos, end=matrices$genomic.pos))
    ## TODO sync chromosome names if not matching between VCFs
    ## GenomeInfoDb::seqlevels(target.sites) <- gsub(GenomeInfoDb::seqlevels(target.sites), pattern = 'chr', replacement = '')
    
    popul.vcf <- popul.vcf.list[[chr]]
    ## Check if VCF file exists along with corresponding TABIX file
    if (file.exists(popul.vcf)) {
      popul.vcf.tbi <- paste0(popul.vcf, '.tbi')
      if (!file.exists(popul.vcf.tbi)) {
        message(paste0("        Inderxing VCF file: ", popul.vcf))
        idx <- Rsamtools::indexTabix(popul.vcf, "vcf")
      }
    } else {
      stop("VCF file: ", popul.vcf, " doesn't exists, quitting ...")
    }
    
    vcf.ranges <- VariantAnnotation::readVcfAsVRanges(x = popul.vcf, param = VariantAnnotation::ScanVcfParam(info = NA, geno = 'GT', which = target.sites))
    genome(vcf.ranges) <- NA
    ## Remove multi nucleotide variants
    vcf.ranges <- vcf.ranges[nchar(VariantAnnotation::ref(vcf.ranges)) == 1 & nchar(VariantAnnotation::alt(vcf.ranges)) == 1]
    ## Remove duplicated positions
    mask <- duplicated(paste0(start(vcf.ranges), '_', VariantAnnotation::sampleNames(vcf.ranges)))
    vcf.ranges <- vcf.ranges[!mask]
    all.1000G.samples <- as.character(unique(VariantAnnotation::sampleNames(vcf.ranges)))
    #vcf.ranges.grl <- split(vcf.ranges, sampleNames(vcf.ranges))
    
    ## Loop over each cell specific WC region and calculate similarity to all 1000G samples
    region.compare <- list()
    for (i in seq_along(matrices$row.IDs)) {
      region.ID <- matrices$row.IDs[i]
      cell.ID <- strsplit(as.character(region.ID), "__")[[1]][1]
      message("    Working on WC region: ", region.ID)
      region.gr <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges::IRanges(start=matrices$genomic.pos, end=matrices$genomic.pos), 
                           crick=matrices$crick.bases[i,], 
                           watson=matrices$watson.bases[i,])
      #GenomeInfoDb::seqlevels(region.gr) <- gsub(GenomeInfoDb::seqlevels(region.gr), pattern = 'chr', replacement = '')
      ## Keep only SNV positions covered either in crick or watson reads
      region.gr <- region.gr[region.gr$crick > 0 | region.gr$watson > 0]
      ## Translate bases
      region.gr$crick <- chartr("01234", "NACGT", region.gr$crick)
      region.gr$watson <- chartr("01234", "NACGT", region.gr$watson)
      
      message("        Comparing SNVs ...")
      
      ## Translate genotypes to bases
      ref <- VariantAnnotation::ref(vcf.ranges)
      alt <- VariantAnnotation::alt(vcf.ranges)
      sampleNames <- VariantAnnotation::sampleNames(vcf.ranges)
      vcf.ranges.gr <- as(vcf.ranges, 'GRanges')
      #alleles <- strsplit(vcf.ranges.gr$GT, "\\/|\\|")
      #alleles <- do.call(rbind, alleles)
      alleles <- stringr::str_split_fixed(string = vcf.ranges.gr$GT, pattern = '\\/|\\|', n = 2)
      allele1 <- ifelse(alleles[,1] == 0, ref, alt)
      allele2 <- ifelse(alleles[,2] == 0, ref, alt)
      sample.gr <- vcf.ranges.gr[,0]
      sample.gr$allele1 <- allele1
      sample.gr$allele2 <- allele2
      sample.gr$sampleNames <- sampleNames 
      ## Keep only positions present in both Strand-seq and 1000G data
      region.gr.sub <- IRanges::subsetByOverlaps(region.gr, sample.gr) 
      sample.gr.sub <- IRanges::subsetByOverlaps(sample.gr, region.gr.sub)
      ## Replicate Strand-seq SNVs by the number of 1000G samples
      region.gr.sub <- rep(region.gr.sub, length(unique(sample.gr.sub$sampleNames)))
      ## Crick comparison
      crick.compare <- data.frame(crick.snv = region.gr.sub$crick,
                                  toAllele1 = region.gr.sub$crick == sample.gr.sub$allele1,
                                  toAllele2 = region.gr.sub$crick == sample.gr.sub$allele2,
                                  sample = sample.gr.sub$sampleNames)
      crick.compare.summary <- crick.compare %>% dplyr::filter(crick.snv != 'N') %>% dplyr::group_by(sample) %>%
        dplyr::summarise(all.pos.crick = dplyr::n(), agree.crick = pmax(length(toAllele1[toAllele1 == TRUE]), length(toAllele2[toAllele2 == TRUE])))
      if (nrow(crick.compare.summary) < length(all.1000G.samples)) {
        crick.compare.summary <- data.frame(sample=all.1000G.samples, all.pos.crick=0, agree.crick=0)
      }
      ## Watson comparison
      watson.compare <- data.frame(watson.snv = region.gr.sub$watson,
                                   toAllele1 = region.gr.sub$watson == sample.gr.sub$allele1,
                                   toAllele2 = region.gr.sub$watson == sample.gr.sub$allele2,
                                   sample = sample.gr.sub$sampleNames)
      watson.compare.summary <- watson.compare %>% dplyr::filter(watson.snv != 'N') %>% dplyr::group_by(sample) %>%
        dplyr::summarise(all.pos.watson = dplyr::n(), agree.watson = pmax(length(toAllele1[toAllele1 == TRUE]), length(toAllele2[toAllele2 == TRUE])))
      if (nrow(watson.compare.summary) < length(all.1000G.samples)) {
        watson.compare.summary <- data.frame(sample=all.1000G.samples, all.pos.watson=0, agree.watson=0)
      }
      
      compare.summary <- cbind(crick.compare.summary, watson.compare.summary[,c(2:3)])
      compare.summary$chr <- chr
      compare.summary$cell.ID <- cell.ID
      region.compare[[i]] <- compare.summary
    }
    compare.per.chr[[length(compare.per.chr) + 1]] <- do.call(rbind, region.compare)
  }
  compare.per.chr.df <- do.call(rbind, compare.per.chr)
  return(compare.per.chr.df)
}
