#' Correct StrandPhaseR phasing over inverted regions
#' 
#' This function attempts to correct phasing over inverted regions using both \pkg{\link{breakpointR}} and \pkg{\link{StrandPhaseR}} results.
#' 
#' @param input.bams A path to a VCF file to be loaded.
#' @param outputfolder A path to a folder where all summary and intermediate results will be stored.
#' @param inv.bed A BED formatted file that contains likely inverted regions to evaluate.
#' @param recall.phased If set to \code{TRUE},inversion signatures will be re-called using phased Strand-seq reads (Useful to redefine HET inversion boundaries).
#' @param het.genotype If set to 'strict' heterozygous inversion has to be supported by both 'breakpointR.data' and 'strandphaseR.data',
#' if set to 'lenient' support from either 'breakpointR.data' only is sufficient.
#' @param chromosomes Limit analysis to a certain chromosomes only.
#' @param snv.positions A path to a VCF files containing SNV positions to phase.
#' @param breakpointR.data A path to results obtained by running \code{\link{breakpointR}} package on a given sample.
#' @param strandphaseR.data A path to results obtained by running \code{\link{StrandPhaseR}} package on a given sample.
#' @param vcfs.files A path to a folder where all VCF files to be corrected are stored (sample specific).
#' @param overwrite.results Set to \code{FALSE} if existing corrected VCF files should be used instead of creating a new copy of original VCF file. 
#' @importFrom breakpointR synchronizeReadDir removeReadPileupSpikes runBreakpointr genotype.binom
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges subsetByOverlaps
#' @importFrom utils read.table write.table
#' @inheritParams bamregion2GRanges
#' @inheritParams phaseChromosome
#' @inheritParams phaseHETinversion
#' @inheritParams vcf2vranges
#' @return \code{NULL} or \code{data.frame} depending if 'outputfolder' is defined.
#' @author David Porubsky
#' @export
#'
correctInvertedRegionPhasing <- function(input.bams, outputfolder=NULL, inv.bed=NULL, recall.phased=FALSE, het.genotype='strict', chromosomes=NULL, snv.positions=NULL, breakpointR.data=NULL, strandphaseR.data=NULL, pairedEndReads=TRUE, min.mapq=10, vcfs.files=NULL, lookup.bp=1000000, lookup.blacklist=NULL, bsGenome=NULL, ref.fasta=NULL, assume.biallelic=TRUE, overwrite.results=TRUE) {
  
  ## Check user input ##
  ######################
  ## Check if set of inverted regions is defined
  if (!is.null(inv.bed) & file.exists(inv.bed)) {
    ## Load set of inverted regions
    inv.df <- utils::read.table(file = inv.bed, header = FALSE, stringsAsFactors = FALSE)
    colnames(inv.df)[c(1:3)] <- c('seqnames', 'start', 'end')
    inv.gr <- makeGRangesFromDataFrame(df = inv.df)
    call.inversion <- FALSE
  } else {
    message("Parameter 'inv.bed' is not specified or given file doesn't exists.\nLikely inverted regions will be called de novo.")
    call.inversion <- TRUE
    recall.phased <- TRUE
  } 
  
  ## Check input format of blacklisted regions and load a BED file if needed
  if (!is.null(lookup.blacklist) & !class(lookup.blacklist) == 'GRanges') {
    if (file.exists(lookup.blacklist)) {
      lookup.blacklist.df <- utils::read.table(file = lookup.blacklist, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
      if (ncol(lookup.blacklist.df) >= 3) {
        lookup.blacklist.df <- lookup.blacklist.df[,c(1:3)]
        colnames(lookup.blacklist.df) <-  c('seqnames', 'start', 'end')
        lookup.blacklist <- GenomicRanges::makeGRangesFromDataFrame(lookup.blacklist.df)
      } else {
        lookup.blacklist <- NULL
        warning("The BED file, '", lookup.blacklist, "' doesn't contain required fields ('chr.name', 'start', 'end').")
      }
    } else {
      lookup.blacklist <- NULL
      warning("The BED file, '", lookup.blacklist, "' doesn't exists.")
    }
  }
  
  ## Abort function execution if breakpointR.data are not defined
  if (is.null(breakpointR.data) & !file.exists(breakpointR.data)) {
    stop("Required parameter 'breakpointR.data' is not defined or file doesn't exists, aborting ...")
  } else {
    ## Load breakpointR.data and create a composite file
    breakp.data <- list.files(breakpointR.data, full.names = TRUE)
    dir.reads <- breakpointR::synchronizeReadDir(files2sync = breakp.data)
  }
  
  ## Abort function execution if strandphaseR.data are not defined
  if (is.null(strandphaseR.data) & !file.exists(strandphaseR.data)) {
    stop("Required parameter 'strandphaseR.data' is not defined or file doesn't exists, aborting ...")
  } else {
    phased.data <- list.files(strandphaseR.data, pattern = "reads.RData", full.names = TRUE)
    phased.files <- basename(basename(phased.data))
    phased.chroms <- sapply(phased.files, function(x) strsplit(x, '_')[[1]][1])
  }
  
  ## Get chromosomes to process
  if (is.null(chromosomes)) {
    chromosomes <- phased.chroms
  }
  chroms2use <- intersect(chromosomes, phased.chroms)
  
  ## Create outputfolder if it doesn't exists
  if (!is.null(outputfolder)) {
    if (!dir.exists(outputfolder)) {
      dir.create(outputfolder)
    }
  }
  
  ## Call inversions in composite file if inverted regions not defined ##
  #######################################################################
  if (call.inversion) {
    ## Remove pileup spikes in coverage
    dir.reads <- breakpointR::removeReadPileupSpikes(dir.reads, max.pileup = 10)
    ## Analyze breakpoints in composite files
    breakpoints <- breakpointR::runBreakpointr(bamfile = dir.reads, 
                                               pairedEndReads = pairedEndReads, 
                                               chromosomes = chromosomes,
                                               windowsize = 5000, 
                                               peakTh = 0.33, 
                                               binMethod = "size", 
                                               genoT = 'binom',
                                               background = 0.05, 
                                               minReads = 50)
    break.regions <- breakpoints$counts
    inv.gr <- break.regions[break.regions$states != 'cc']
  }
  
  ## Process each chromosome separately ##
  ########################################
  correction.summary <- list()
  for (chr in chroms2use) {
    ## Initialize corrected VCF file
    if (!is.null(vcfs.files)) {
      vcf.filename <- paste0(chr, '_phased.vcf')
      vcf.file <- file.path(vcfs.files, vcf.filename)
      if (file.exists(vcf.file)) {
        ## Make a copy of an input VCF file
        vcf.file.corr <- file.path(dirname(vcf.file), gsub(vcf.filename, pattern = '\\.vcf', replacement = '_INVcorr.vcf', ignore.case = TRUE))
        if (!file.exists(vcf.file.corr) | overwrite.results == TRUE) {
          message("Writing a VCF file copy, '", vcf.file.corr,"' for inversion correction.")
          file.copy(vcf.file, vcf.file.corr, overwrite = TRUE)
        } else {
          message("Using existing VCF file, '", vcf.file.corr,"' for inversion correction.")
        } 
      } else {
        warning(paste0("VCF file: ", vcf.file, " doesn't exists, skipping correction ..."))
      }  
    } else {
      warning("Parameter 'vcfs.files' not defined, skipping correction ...")
    }
    
    ## Load strandphaseR.data and create a composite phased file (per chromosome)
    phased.data.chr <- file.path(strandphaseR.data, paste0(chr, '_reads.RData'))
    reads <- get(load(phased.data.chr))
    strand(reads$hap1) <- "+" ## HAP1 == CRICK
    strand(reads$hap2) <- "-" ## HAP2 == WATSON
    phased.reads <- c(reads$hap1, reads$hap2)
    phased.reads <- GenomicRanges::sort(phased.reads, ignore.strand=TRUE)
    
    if (recall.phased) {
      breakp.phased <- breakpointR::runBreakpointr(bamfile = phased.reads, 
                                                   pairedEndReads = FALSE, 
                                                   chromosomes = chr,
                                                   windowsize = 10000, 
                                                   peakTh = 0.33, 
                                                   binMethod = "size", 
                                                   genoT = 'binom',
                                                   background = 0.05, 
                                                   minReads = 50)
      detected.ranges <- breakp.phased$counts
      ## Report regions that are genotyped as purely WW or CC as those point to likely HET inversions
      het.inv.gr <- detected.ranges[detected.ranges$states %in% c('ww', 'cc')]
    }
    
    ## Subset inversions for a given chromosome
    inv.gr.chr <- inv.gr[GenomeInfoDb::seqnames(inv.gr) == chr] 
    if (length(inv.gr.chr) > 0) {
      message("Correcting phasing for inversions on chromosome: ", chr)
    } else {
      message("No inversions reported for chromosome: ", chr)
      next
    }
    
    ## Process each inversion separately ##
    #######################################
    for (j in seq_along(inv.gr.chr)) {
      roi.gr <- inv.gr.chr[j]
      ## Get directional and phased reads for a given region
      roi.dir.reads <- IRanges::subsetByOverlaps(dir.reads, roi.gr)
      ## If 'recall.phased' set to TRUE try to use ranges defined by breakpoint calling on phased reads
      if (recall.phased) {
        gr <- IRanges::subsetByOverlaps(het.inv.gr, roi.gr) 
        if (length(gr) > 0) {
          roi.phased.reads <- IRanges::subsetByOverlaps(phased.reads, range(gr))
        } else {
          roi.phased.reads <- IRanges::subsetByOverlaps(phased.reads, roi.gr)
        }
      } else {
        roi.phased.reads <- IRanges::subsetByOverlaps(phased.reads, roi.gr)
      }  
      ## Genotype region of interest (ROI) ##
      ## Remove pileup spikes in coverage
      roi.dir.reads <- breakpointR::removeReadPileupSpikes(gr = roi.dir.reads, max.pileup = 10)
      ## Genotype roi directional reads
      wReads <- roi.dir.reads[GenomicRanges::strand(roi.dir.reads) == '-']
      cReads <- roi.dir.reads[GenomicRanges::strand(roi.dir.reads) == '+']
      dir.reads.genot <- breakpointR::genotype.binom(wReads=length(wReads), cReads=length(cReads), background=0.05, minReads=10, log=TRUE)
      ## Genotype roi phased reads
      wReads <- roi.phased.reads[strand(roi.phased.reads) == '-']
      cReads <- roi.phased.reads[strand(roi.phased.reads) == '+']
      phased.reads.genot <- breakpointR::genotype.binom(wReads=length(wReads), cReads=length(cReads), background=0.05, minReads=10, log=TRUE)
      ## Get inversion genotype
      if (!is.na(dir.reads.genot$bestFit) & !is.na(phased.reads.genot$bestFit)) {
        ## Defined if inversion is homozygous (HOM)
        if (dir.reads.genot$bestFit == 'ww' & phased.reads.genot$bestFit == 'wc') {
          inv.genot <- 'HOM'
        } else if (dir.reads.genot$bestFit == 'wc' & phased.reads.genot$bestFit == 'wc') {
          inv.genot <- 'complex'
        } else {
          inv.genot <- 'REF'
        }
        ## Defined if inversion is heterozygous (HET)
        if (het.genotype == 'strict') {
          if (dir.reads.genot$bestFit == 'wc' & phased.reads.genot$bestFit %in% c('ww', 'cc')) {
            inv.genot <- 'HET'
          }
        } else if (het.genotype == 'lenient') {
          if (dir.reads.genot$bestFit == 'wc') {
            inv.genot <- 'HET'
          }
        }  
      } else {
        inv.genot <- 'lowReads'
      }  
      
      ## Phase HET inversion
      if (inv.genot == 'HET') {
        ## Set range to be phased
        if (recall.phased) {
          het.gr <- IRanges::subsetByOverlaps(het.inv.gr, roi.gr) 
          if (length(gr) > 0) {
            het.gr <- range(het.gr)
          } else {
            het.gr <- roi.gr
          }
        } else {
          het.gr <- roi.gr
        }
        
        het.haps <- phaseHETinversion(input.bams = input.bams, snv.positions=snv.positions, phase.gr=het.gr, lookup.bp = lookup.bp, lookup.blacklist = lookup.blacklist, pairedEndReads = pairedEndReads, min.mapq = min.mapq, bsGenome = bsGenome, ref.fasta = ref.fasta, assume.biallelic = assume.biallelic)
        ## Assign inverted haplotype to either H1 or H2 based on 'phased.reads.genot'
        if (phased.reads.genot$bestFit == 'cc') {
          if (!is.null(het.haps)) {
            het.haps$H1 <- het.haps$ref.phase
            het.haps$H2 <- het.haps$inv.phase
          }
          H1 <- 'ref'
          H2 <- 'inv'
        } else if (phased.reads.genot$bestFit == 'ww') {
          if (!is.null(het.haps)) {
            het.haps$H1 <- het.haps$inv.phase
            het.haps$H2 <- het.haps$ref.phase
          }  
          H1 <- 'inv'
          H2 <- 'ref'
        } else {
          ## Not able to determine keep H1 as REF and H2 as ALT
          if (!is.null(het.haps)) {
            het.haps$H1 <- het.haps$ref.phase
            het.haps$H2 <- het.haps$inv.phase
          }  
          H1 <- '.'
          H2 <- '.'
        }
      }
      
      ## Correct haplotypes overlapping with inverted region (roi.gr)
      #if (!is.null(vcfs.files)) {
      #vcf.filename <- paste0(chr, '_phased.vcf')
      #vcf.file <- file.path(vcfs.files, vcf.filename)
      if (file.exists(vcf.file.corr)) {
        if (inv.genot == 'HOM') {
          correctHomInv(correct.gr = roi.gr, vcf.file = vcf.file.corr)
          H1 <- 'inv'
          H2 <- 'inv'
        } else if (inv.genot == 'HET') {
          if (!is.null(het.haps)) {
            correctHetInv(correct.gr = het.gr, vcf.file = vcf.file.corr, het.haps = het.haps)
          }
        } else {
          H1 <- inv.genot
          H2 <- inv.genot
        }
      } else {
        warning(paste0("VCF file: ", vcf.file, " doesn't exists, skipping correction ..."))
      }  
      #} else {
      #  warning("Parameter 'vcfs.files' not defined, skipping correction ...")
      #}
      ## Report INV correction summary
      if (inv.genot == 'HET') { 
        het.region <- as.character(het.gr)
      } else {
        het.region <- 'NA'
      }
      report.df <- data.frame(bed.region=as.character(roi.gr), 
                              dir.reads.gt=dir.reads.genot$bestFit, 
                              phased.reads.gt=phased.reads.genot$bestFit, 
                              predict.inv.gt=inv.genot,
                              het.region=het.region, 
                              H1=H1, 
                              H2=H2)
      correction.summary[[length(correction.summary) + 1]] <- report.df
    } ## End of INV loop
  } ## End of CHR loop
  ## Report correction summary table
  correction.summary.df <- do.call(rbind, correction.summary)
  if (!is.null(outputfolder)) {
    destination <- file.path(outputfolder, 'invPhasingCorrection.summary.tsv')
    utils::write.table(correction.summary.df, file = destination, quote = FALSE, sep = '\t', row.names = FALSE)
    return(NULL)
  } else {
    return(correction.summary.df)
  }  
} ## End of function  


#' Phase a heterozygous inversion using Strand-seq data
#' 
#' This function takes as an input region deemed to be a heterozygous inversion and attempts to phase SNVs inside this region using
#' haplotype informative Strand-seq reads from multiple single cells.
#' 
#' @param phase.gr A \code{\link{GRanges}} object of region to be phased. This region should point to a heterozygous inversion site.
#' @param lookup.bp A number of nucleotides, downstream and upstream, from the heterozygous inversion site ('phase.gr') to be genotyped.
#' @param lookup.blacklist A \code{\link{GRanges}} object or a path to a BED file containing a set of ranges to be excluded 
#' when extending 'phase.gr' by 'lookup.bp'. The total size of 'lookup.bp' is kept after filtering.
#' @param assume.biallelic If set to \code{TRUE} parameter 'snv.positions' is expected to contain biallelic loci (0/1, 1/0) and thus
#' gaps in haplotypes will be filled accordingly.
#' @param verbose Is set to \code{TRUE} function will provide a more detailed messaging of ongoing analysis steps.
#' @importFrom bamsignals bamCount
#' @importFrom Biostrings alphabetFrequency Views
#' @importFrom IRanges findOverlaps
#' @importFrom breakpointR genotype.binom
#' @importFrom GenomeInfoDb keepSeqlevels
#' @inheritParams correctInvertedRegionPhasing
#' @inheritParams bamregion2GRanges
#' @inheritParams phaseChromosome
#' @inheritParams vcf2vranges
#' @inheritParams exportVCF
#' @return A \code{data.frame} with phased alleles per inverted and reference haplotype.
#' @author David Porubsky
#' @export
#'
phaseHETinversion <- function(input.bams=NULL, snv.positions=NULL, phase.gr=NULL, lookup.bp=1000000, lookup.blacklist=NULL, pairedEndReads=TRUE, min.mapq=10, bsGenome=NULL, ref.fasta=NULL, assume.biallelic=TRUE, verbose=FALSE) {
  
  ## Load BSgenome
  if (class(bsGenome) != 'BSgenome') {
    if (is.character(bsGenome)) {
      suppressPackageStartupMessages(library(bsGenome, character.only=TRUE))
      bsGenome <- eval(parse(text=bsGenome)) # replacing string by object
    } else {
      bsGenome <- NULL
    }
  }
  
  ## Check if ref.fasta file is in a valid FASTA format
  if (!is.null(ref.fasta)) {
    if (class(Rsamtools::FaFile(ref.fasta)) != "FaFile") {
      ref.fasta <- NULL
      warning("User defined reference FASTA file, '", ref.fasta, "' is not in proper FASTA format!!!")
    }
  }
  
  ## Get chromosome name from 'phase.gr' object
  chromosome <- as.character(unique(seqnames(phase.gr)))
  
  ## Get SNV positions for a region defined in 'phase.gr'
  if (!is.null(snv.positions)) {
    snvs <- suppressMessages( vcf2ranges(vcfFile=snv.positions, genotypeField=1, chromosome = chromosome) )
    region.snvs <- IRanges::subsetByOverlaps(snvs, phase.gr)
    ## Remove duplicated SNV positions
    region.snvs <- region.snvs[!duplicated(region.snvs)]
  } else {
    stop("Required parameter 'snv.positions' not defined!!!")
  }
  ## Report empty data.frame object if there are no SNVs overlapping region of interest
  if (length(region.snvs) > 0) {
    
    ## Loop over all BAM files and select cells with WC genotype over region of interest
    file.list <- list.files(path = input.bams, pattern = '\\.bam$', full.names = TRUE) 
    
    inv.reads <- GenomicRanges::GRangesList()
    ref.reads <- GenomicRanges::GRangesList()
    for (i in seq_along(file.list)) {
      bamfile <- file.list[i]
      bamfile.name <- basename(bamfile)
      if (verbose) {
        message("Processing bamfile: ", bamfile.name) 
      }  
      ## Genotype inversions flanks ##
      lookup.gr <- GenomicRanges::resize(phase.gr, width = (2 * lookup.bp) + width(phase.gr), fix = 'center')
      ## Make sure lookup.gr is no larger than set seqlength
      lookup.gr <- GenomicRanges::trim(lookup.gr)
      flanks.gr <- GenomicRanges::setdiff(lookup.gr, phase.gr)
      ## Account for blacklisted regions
      if (!is.null(lookup.blacklist)) {
        downANDupsteam <- GenomicRanges::gaps(phase.gr)
        downANDupsteam <- downANDupsteam[strand(downANDupsteam ) == '*']
        unique <- suppressWarnings( GenomicRanges::setdiff(downANDupsteam, lookup.blacklist) )
        hits <- IRanges::findOverlaps(unique, downANDupsteam)
        unique.grl <- GenomicRanges::split(unique, subjectHits(hits))
        flanks.grl <- GenomicRanges::GRangesList()
        for (j in seq_along(unique.grl)) {
          gr <- unique.grl[[j]]
          if (j == 1) {
            flanks.grl[[j]] <- sliceRanges(gr, slice.size = lookup.bp, from = 'end')
          } else {
            flanks.grl[[j]] <- sliceRanges(gr, slice.size = lookup.bp, from = 'start')
          }
        }
        flanks.gr <- unlist(flanks.grl)
      }
      ## Set parameter for bamsignals counts
      paired.end <- 'ignore'
      if (pairedEndReads) {
        paired.end <- 'filter'
      }
      max.frag <- 1000
      counts <- bamsignals::bamCount(bamfile, flanks.gr, mapq=min.mapq, filteredFlag=1024, paired.end=paired.end, tlenFilter=c(0, max.frag), verbose=FALSE, ss=TRUE)
      flanks.genot <- breakpointR::genotype.binom(wReads=sum(counts['antisense',]), cReads=sum(counts['sense',]), background=0.05, minReads=50, log=TRUE)
      
      if (flanks.genot$bestFit %in% c('ww','cc')) {
        ## Load reads from inverted region
        bam.reads <- bamregion2GRanges(bamfile, region=phase.gr, pairedEndReads=pairedEndReads, min.mapq=min.mapq)
      } else {
        if (verbose) {
          message("    Uninformative cell, skipping ...")
        }  
        next
      } 
      
      ## Store phased reads
      if (flanks.genot$bestFit == 'ww') {
        inv.reads[[length(inv.reads) + 1]] <- bam.reads[strand(bam.reads) == '+']
        ref.reads[[length(ref.reads) + 1]] <- bam.reads[strand(bam.reads) == '-']
      } else if (flanks.genot$bestFit == 'cc') {
        inv.reads[[length(inv.reads) + 1]] <- bam.reads[strand(bam.reads) == '-']
        ref.reads[[length(ref.reads) + 1]] <- bam.reads[strand(bam.reads) == '+']
      }
    }
    inv.reads <- unlist(inv.reads, use.names = FALSE)
    ref.reads <- unlist(ref.reads, use.names = FALSE)
    
    if (length(inv.reads) != 0 & length(ref.reads) != 0) {
      ## Make sure only reads from a chromosome in phase.gr are kept
      inv.reads <- GenomeInfoDb::keepSeqlevels(inv.reads, value = chromosome, pruning.mode = 'coarse')
      ref.reads <- GenomeInfoDb::keepSeqlevels(ref.reads, value = chromosome, pruning.mode = 'coarse')
      
      ## Extract alleles at variable positions ##
      ## Extract read sequences for watson and crick reads
      inv.reads.seq <- GenomicRanges::mcols(inv.reads)$seq
      ref.reads.seq <- GenomicRanges::mcols(ref.reads)$seq
      ## Extract base qualities for watson and crick reads
      #inv.reads.qual <- mcols(inv.reads)$qual
      #ref.reads.qual <- mcols(ref.reads)$qual
      
      ## get piles of bases at each variable position
      piles.inv.reads <- GenomicAlignments::pileLettersAt(inv.reads.seq, seqnames(inv.reads), start(inv.reads), mcols(inv.reads)$cigar, region.snvs)
      piles.ref.reads <- GenomicAlignments::pileLettersAt(ref.reads.seq, seqnames(ref.reads), start(ref.reads), mcols(ref.reads)$cigar, region.snvs)
      
      ## get piles of qualities at each variable position
      #quals.inv.reads <- GenomicAlignments::pileLettersAt(inv.reads.qual, seqnames(inv.reads), start(inv.reads), mcols(inv.reads)$cigar, region.snvs)
      #quals.ref.reads <- GenomicAlignments::pileLettersAt(ref.reads.qual, seqnames(ref.reads), start(ref.reads), mcols(ref.reads)$cigar, region.snvs)
      
      if (sum(width(piles.inv.reads)) > 0) {
        ## Get inverted alleles indices
        inv.bases.freq <- Biostrings::alphabetFrequency(piles.inv.reads, baseOnly=TRUE)
        ## Remove empty sites
        inv.bases.idx <- which(rowSums(inv.bases.freq) > 0)
        inv.bases.freq <- inv.bases.freq[drop=FALSE, inv.bases.idx,]
        # Remove positions with ambiguous allele counts (max and second max equally abundant)
        filt.ambig <- apply(inv.bases.freq, 1, function(x) (sum(x) - max(x)) < max(x))
        inv.bases.freq <- inv.bases.freq[drop=FALSE, filt.ambig,]
        inv.bases.idx <-  inv.bases.idx[filt.ambig]
        # Get max covered alleles
        max.idx <- apply(inv.bases.freq, 1, which.max)
        # Keep only standard nucleotides
        filt.nonstand <- max.idx <= 4
        max.idx <- max.idx[filt.nonstand]
        inv.bases.idx <- inv.bases.idx[filt.nonstand]
        inv.bases <- chartr("1234", "ACGT", max.idx)
      } else {
        inv.bases <- 'N'
        inv.bases.idx <- 1
      }  
      
      # if (sum(width(piles.ref.reads)) > 0) {
      #   ## Get reference alleles indices
      #   ref.bases.freq <- Biostrings::alphabetFrequency(piles.ref.reads, baseOnly=TRUE)
      #   ref.bases.idx <- which( ref.bases.freq > 0, arr.ind=TRUE)
      #   ref.bases.idx <- ref.bases.idx[drop=FALSE, order(ref.bases.idx[,1]),] # sort by SNV position
      #   ref.bases.idx <- ref.bases.idx[drop=FALSE, ref.bases.idx[, 2] <= 4, ] # keep only standard nucleotides
      #   ref.bases <- chartr("1234", "ACGT", ref.bases.idx[, 2])
      # } else {
      #   ref.bases <- 'N'
      #   ref.bases.idx <- matrix(c(1,1), nrow = 1)
      # }  
      
      if (sum(width(piles.ref.reads)) > 0) {
        ## Get inverted alleles indices
        ref.bases.freq <- Biostrings::alphabetFrequency(piles.ref.reads, baseOnly=TRUE)
        ## Remove empty sites
        ref.bases.idx <- which(rowSums(ref.bases.freq) > 0)
        ref.bases.freq <- ref.bases.freq[drop=FALSE, ref.bases.idx,]
        # Remove positions with ambiguous allele counts (max and second max equally abundant)
        filt.ambig <- apply(ref.bases.freq, 1, function(x) (sum(x) - max(x)) < max(x))
        ref.bases.freq <- ref.bases.freq[drop=FALSE, filt.ambig,]
        ref.bases.idx <-  ref.bases.idx[filt.ambig]
        # Get max covered alleles
        max.idx <- apply(ref.bases.freq, 1, which.max)
        # Keep only standard nucleotides
        filt.nonstand <- max.idx <= 4
        max.idx <- max.idx[filt.nonstand]
        ref.bases.idx <- ref.bases.idx[filt.nonstand]
        ref.bases <- chartr("1234", "ACGT", max.idx)
      } else {
        ref.bases <- 'N'
        ref.bases.idx <- 1
      }
      
      ## Assemble haplotypes
      # pos.idx <- unique(sort(c(inv.bases.idx[, 1], ref.bases.idx[, 1])))
      # pos.gen <- start(region.snvs)[pos.idx]
      # inv.allele <- rep('.', times=length(pos.gen))
      # inv.allele[match(inv.bases.idx[, 1], pos.idx)] <- inv.bases
      # ref.allele <- rep('.', times=length(pos.gen))
      # ref.allele[match(ref.bases.idx[, 1], pos.idx)] <- ref.bases
      
      pos.idx <- unique(sort(c(inv.bases.idx, ref.bases.idx)))
      pos.gen <- start(region.snvs)[pos.idx]
      inv.allele <- rep('.', times=length(pos.gen))
      inv.allele[match(inv.bases.idx, pos.idx)] <- inv.bases
      ref.allele <- rep('.', times=length(pos.gen))
      ref.allele[match(ref.bases.idx, pos.idx)] <- ref.bases
      
      ## Assign reference and alternative alleles
      if (!is.null(bsGenome)) {
        ## Extract reference alleles from the reference bsgenome object if available
        snv.ranges <- GenomicRanges::GRanges(seqnames=chromosome, IRanges(start=pos.gen, end=pos.gen))
        ref.base <- Biostrings::Views(bsGenome, snv.ranges)
        ref.base <- as(ref.base, "DNAStringSet")
        ref.base <- as(ref.base, "vector")
      } else if (!is.null(ref.fasta)) {
        ## Extract SNV bases from user defined reference FASTA file
        fa.file <- open(Rsamtools::FaFile(ref.fasta))
        snv.seq <- Rsamtools::scanFa(file = fa.file, param = snv.ranges, as = "DNAStringSet")     
        names(snv.seq) <- NULL
        ref.base <- as.character(snv.seq)
      } else {
        ref.base <- rep('N', length(pos.gen))
      }
        
      ## Report inverted and reference haplotypes
      inv.phase <- rep('.', times=length(pos.gen))
      ref.phase <- rep('.', times=length(pos.gen))
      inv.phase[inv.allele == ref.base] <- 0
      inv.phase[inv.allele != ref.base & inv.allele != '.'] <- 1
      ref.phase[ref.allele == ref.base] <- 0
      ref.phase[ref.allele != ref.base & ref.allele != '.'] <- 1
      inv.allele[inv.allele == '.'] <- 'N'
      ref.allele[ref.allele == '.'] <- 'N'
      
      if (assume.biallelic) {
        inv.phase.gaps <- which(inv.phase == '.')
        inv.phase[inv.phase.gaps] <- dplyr::recode(ref.phase[inv.phase.gaps], '0' = '1', '1' = '0')
        ref.phase.gaps <- which(ref.phase == '.')
        ref.phase[ref.phase.gaps] <- dplyr::recode(inv.phase[ref.phase.gaps], '0' = '1', '1' = '0')
        inv.allele[inv.phase == '0'] <- ref.base[inv.phase == '0']
        ref.allele[ref.phase == '0'] <- ref.base[ref.phase == '0']
      }  
      phased.df <- data.frame(pos.gen=pos.gen, ref.allele=ref.allele, inv.allele=inv.allele, ref.phase=ref.phase, inv.phase=inv.phase, ref.base=ref.base)
    } else {
      phased.df <- NULL
    }
  } else {
    phased.df <- NULL
  }
  return(phased.df)
}


#' Correct homozygous inversion in StrandPhaseR VCF file.
#' 
#' This function takes as an input region deemed to be a homozygous inversion and flips H1 and H2 haplotypes for SNVs that 
#' overlap with this region.
#' 
#' @param correct.gr A \code{\link{GRanges}} object of homozygous inversion site to be corrected.
#' @param vcf.file A \pkg{\link{StrandPhaseR}} formatted VCF file to be corrected for inversion phasing.
#' @param ID A unique id to be appended at the end of each corrected VCF file.
#' @importFrom VariantAnnotation readVcfAsVRanges ScanVcfParam
#' @inheritParams correctInvertedRegionPhasing
#' @return \code{NULL}
#' @author David Porubsky
#' @export
#'
correctHomInv <- function(correct.gr=NULL, vcf.file=NULL, ID='') {
  ## Read VCF file header
  file.conn <- file(vcf.file, "r")
  line <- readLines(file.conn, 1)
  header <- list()
  while(grepl("#", line)) {
    header[[length(header ) + 1]] <- line
    line <- readLines(file.conn, 1)
  }
  #header <- cat(unlist(header), sep = '\n')
  
  ## Read VCF file
  #vcf.data <- read.table(vcf.file, header = FALSE, comment.char = "#")
  suppressWarnings( vcf.data <- VariantAnnotation::readVcfAsVRanges(x = vcf.file, 
                                                                    param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA, geno = c('GT','Q1','Q2','P1','P2') )) )
  ## Filter positions with no reference allele and two alternative alleles
  vcf.data <- vcf.data[!vcf.data$GT %in% c('1|2', '1/2')]
  ## Make sure that missing values are converted from NAs to dots
  vcf.data$Q1[is.na(vcf.data$Q1)] <- '0'
  vcf.data$Q2[is.na(vcf.data$Q2)] <- '0'
  vcf.data$P1[is.na(vcf.data$P1)] <- '0'
  vcf.data$P2[is.na(vcf.data$P2)] <- '0'
  
  ## Flip haplotypes overlapping HOM inversion
  hits <- IRanges::findOverlaps(vcf.data, correct.gr, type = 'within')
  snv.idx <- S4Vectors::queryHits(hits)
  gt.toFlip <- vcf.data$GT[snv.idx]
  gt.toFlip[vcf.data$GT[snv.idx] == '1|0'] <- '0|1'
  gt.toFlip[vcf.data$GT[snv.idx] == '0|1'] <- '1|0'
  vcf.data$GT[snv.idx] <- gt.toFlip
  
  ## Print VCF ##
  if (is.character(ID) & nchar(ID) > 0) {
    destination <- gsub(vcf.file, pattern = '\\.vcf', replacement = paste0('_', ID, '\\.vcf'))
  } else {
    destination <- vcf.file
  }  
  ## Print header
  cat(unlist(header), sep = '\n', file = destination)
  ## Print corrected phasing
  gt <- apply(mcols(vcf.data), 1, function(x) paste(x, collapse = ':'))
  print.df <- data.frame(chr=unique(seqnames(vcf.data)),
                         pos=start(vcf.data),
                         id=rep('.', length(vcf.data)),
                         ref=VariantAnnotation::ref(vcf.data),
                         alt=VariantAnnotation::alt(vcf.data),
                         qual=rep('.', length(vcf.data)),
                         filter=rep('.', length(vcf.data)),
                         info=rep('.', length(vcf.data)),
                         format='GT:Q1:Q2:P1:P2',
                         gt=gt)
  utils::write.table(print.df, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
  return(NULL)
}


#' Correct heterozygous inversion in StrandPhaseR VCF file.
#' 
#' This function takes as an input region deemed to be a heterozygous inversion along with phased SNVs within a region and uses
#' this information to add/correct into the \pkg{\link{StrandPhaseR}} formatted VCF file.
#' 
#' @param het.haps A \code{data.frame} with phased alleles per inverted and reference haplotype.
#' @importFrom VariantAnnotation readVcfAsVRanges ScanVcfParam
#' @inheritParams correctInvertedRegionPhasing
#' @inheritParams correctHomInv
#' @return \code{NULL}
#' @author David Porubsky
#' @export
#'
correctHetInv <- function(correct.gr=NULL, vcf.file=NULL, het.haps=NULL, ID='') {
  ## Check if het.haps parameter is defined and that it contains phase information
  if (!is.null(het.haps) & any(grepl(colnames(het.haps), pattern = 'phase'))) {
    ## Read VCF file header
    file.conn <- file(vcf.file, "r")
    line <- readLines(file.conn, 1)
    header <- list()
    while(grepl("#", line)) {
      header[[length(header ) + 1]] <- line
      line <- readLines(file.conn, 1)
    }
    #header <- cat(unlist(header), sep = '\n')
    
    ## Read VCF file
    suppressWarnings( vcf.data <- VariantAnnotation::readVcfAsVRanges(x = vcf.file, 
                                                                      param = VariantAnnotation::ScanVcfParam(fixed=c('ALT'), info = NA, geno = c('GT','Q1','Q2','P1','P2') )) )
    ## Filter positions with no reference allele and two alternative alleles
    vcf.data <- vcf.data[!vcf.data$GT %in% c('1|2', '1/2')]
    ## Make sure that missing values are converted from NAs to zeroes
    vcf.data$Q1[is.na(vcf.data$Q1)] <- '0'
    vcf.data$Q2[is.na(vcf.data$Q2)] <- '0'
    vcf.data$P1[is.na(vcf.data$P1)] <- '0'
    vcf.data$P2[is.na(vcf.data$P2)] <- '0'
    
    ## Remove SNVs overlapping with HET inversion
    hits <- IRanges::findOverlaps(vcf.data, correct.gr, type = 'within')
    snv.idx <- S4Vectors::queryHits(hits)
    if (length(snv.idx) > 0) {
      vcf.data <- vcf.data[-snv.idx]
    }
    
    #gt <- paste0(het.haps$ref.phase, '|', het.haps$inv.phase) ## Wrongly assuming H1 is always a ref.phase ... !!!
    gt <- paste0(het.haps$H1, '|', het.haps$H2)
    gt <- paste0(gt, ':0:0:0:0')
    ## Get REF and ALT alleles
    REF <- het.haps$ref.base
    ALT <- rep('N', length(gt))
    ALT[het.haps$ref.phase == 1] <- het.haps$ref.allele[het.haps$ref.phase == 1]
    ALT[het.haps$inv.phase == 1] <- het.haps$inv.allele[het.haps$inv.phase == 1]  
    add.df <- data.frame(chr=unique(seqnames(vcf.data)),
                         pos=het.haps$pos.gen,
                         id=rep('.', nrow(het.haps)),
                         ref=REF,
                         alt=ALT,
                         qual=rep('.', nrow(het.haps)),
                         filter=rep('.', nrow(het.haps)),
                         info=rep('.', nrow(het.haps)),
                         format='GT:Q1:Q2:P1:P2',
                         gt=gt)
    
    ## Print VCF ##
    if (is.character(ID) & nchar(ID) > 0) {
      destination <- gsub(vcf.file, pattern = '\\.vcf', replacement = paste0('_', ID, '\\.vcf'))
    } else {
      destination <- vcf.file
    }  
    ## Print header
    cat(unlist(header), sep = '\n', file = destination)
    ## Print corrected phasing
    gt <- apply(mcols(vcf.data), 1, function(x) paste(x, collapse = ':'))
    print.df <- data.frame(chr=unique(seqnames(vcf.data)),
                           pos=start(vcf.data),
                           id=rep('.', length(vcf.data)),
                           ref=VariantAnnotation::ref(vcf.data),
                           alt=VariantAnnotation::alt(vcf.data),
                           qual=rep('.', length(vcf.data)),
                           filter=rep('.', length(vcf.data)),
                           info=rep('.', length(vcf.data)),
                           format='GT:Q1:Q2:P1:P2',
                           gt=gt)
    ## Add newly phased variants from HET inversion
    print.df <- rbind(print.df, add.df)
    print.df <- print.df[order(print.df$pos),]
    utils::write.table(print.df, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
  } else {
    stop("Parameter 'het.haps' is not defined or it doesn't contain phase information expected to be in column 4 and 5 !!!")
  }  
}  

## Outside functions ##
# genotype.binom <- function(wReads, cReads, background=0.05, minReads=10, log=FALSE) {
#   ## Set parameters
#   roiReads <- wReads + cReads
#   alpha <- background
#   ## Calculate binomial probabilities for given read counts
#   result <- list(bestFit=NA, pval=NA)
#   if (length(roiReads)==0) {
#     return(result)
#   }
#   if (is.na(roiReads)) {
#     return(result)
#   }
#   if ( roiReads >= minReads ) {
#     ## Test if a given read counts are WW
#     WWpVal <- stats::dbinom(wReads, size = roiReads, prob = 1-alpha, log = log)
#     ## Test if a given read counts are CC
#     CCpVal <- stats::dbinom(wReads, size = roiReads, prob = alpha, log = log)
#     ## Test if a given read counts are WC
#     WCpVal <- stats::dbinom(wReads, size = roiReads, prob = 0.5, log = log)
#     ## Export results
#     pVal <- c(wc = WCpVal, cc = CCpVal, ww = WWpVal)
#     result <- list(bestFit = names(pVal)[which.max(pVal)], pval = max(pVal))
#     return(result)
#   } else {
#     return(result)
#   }
# }


#' Slice a certain length from a set of genomic ranges. 
#' 
#' This function takes as an input a set of ranges stored in \code{\link{GRanges}} object and export slice of certain length
#' starting either from the first or last range.
#' 
#' @param gr A \code{\link{GRanges}} object containing a set of genomic ranges.
#' @param slice.size A total length of sliced ranges to be reported given the starting range defined by 'from'.
#' @param from Set to 'start' or 'end' if to start counting bases to slice from the first or the last range.
#' @return A \code{\link{GRanges}} object.
#' @author David Porubsky
#' @export
#'
sliceRanges <- function(gr, slice.size=1000000, from='start') {
  ## Make sure submitted set of ranges is sorted by position
  gr <- GenomicRanges::sort(gr)
  if (sum(width(gr)) > slice.size) {
    if (from == 'start') {
      ## Get cumulative sum of ranges width
      cs.width <- cumsum(as.numeric(width(gr)))
      cut.idx <- findInterval(cs.width, slice.size)
      cut.idx <- which(cut.idx == 1)[1]
      gr.sub <- gr[1:cut.idx]
      ## Get cut position
      pos <- mapply(seq, start(gr.sub), end(gr.sub))
      pos <- unlist(pos)[1:slice.size]
      cut.pos <- pos[length(pos)]
      ## Adjust size of the last range
      end(gr.sub[length(gr.sub)]) <- cut.pos
    } else if (from == 'end') {
      ## Get cumulative sum of ranges width
      gr.rev <- rev(gr)
      cs.width <- cumsum(as.numeric(width(gr.rev)))
      cut.idx <- findInterval(cs.width, slice.size)
      cut.idx <- which(cut.idx == 1)[1]
      gr.sub <- gr.rev[1:cut.idx]
      ## Get cut position
      pos <- mapply(seq, end(gr.sub), start(gr.sub))
      pos <- unlist(pos)[1:slice.size]
      cut.pos <- pos[length(pos)]
      ## Adjust size of the last range
      start(gr.sub[length(gr.sub)]) <- cut.pos
    } else {
      stop("Set parameter 'from' to either 'start' or 'end'!!!")
    }
  } else {
    gr.sub <- gr
  }  
  return(gr.sub)
}
