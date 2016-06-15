#' Wrapper function for StrandPhaseR
#'
#' This function will move through .bam files in a folder and perform several steps (see Details).
#'
#' 1. extract variable position in WC regions
#' 2. Fill two matrices separately for SNVs found in Watson and Crick reads
#' 3. Sort matrices in order each column in each matrix has lowest amount of conflicting bases
#' 4. Exclude rows/cells which cannot be reliably assigned to only one matrix consensus 
#' 5. For successfully phased rows/cell export W and C reads as a separate haplotype specifiv GRanges object
#'
#' @param bamfilespath Path to the bam files to process
#' @param dataDirectory Output directory. If non-existent it will be created.
#' @param positions Filename with listed position of SNVs for given chromosome (format: chrName SNVpos).
#' @param WCregions Filename of all WC region for a given chromosome (format: chrName:Start:End:FileName).
#' @param chromosomes If only a subset of the chromosomes should be processed, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param min.baseq Minimum base quality to consider a base for phasing.
#' @param num.iterations Number of iteration to sort watson and crick matrices.
#' @param translateBases
#' @param score2qual
#' @param fillGaps 
#' @param splitPhasedReads Set to \code{TRUE} if you want to split reads per haplotype.
#' @param callBreaks
#' @param exportVCF
#' @param bsGenome
 
#' @author David Porubsky
#' @export

runStrandPhaseR <- function(bamfilespath, dataDirectory='./StrandPhaseR_analysis', positions=NULL, WCregions=NULL, chromosomes=NULL, pairedEndReads=TRUE, min.mapq=10, min.baseq=30, num.iterations=2, translateBases=TRUE, score2qual=FALSE, fillMissAllele=NULL, splitPhasedReads=FALSE, callBreaks=FALSE, exportVCF=NULL, bsGenome=NULL) {
  
  ## Check user input
  if (!is.null(exportVCF) & is.null(bsGenome)) {
	warning("VCf file cannot be created because reference genome is NULL")
  }	
	
  ## Creating directories for data export
  if (!file.exists(dataDirectory)) {
    dir.create(dataDirectory)
  }
  
  phased.store <- file.path(dataDirectory, 'Phased')
  if (!file.exists(phased.store)) {
    dir.create(phased.store)
  }
  
  data.store <- file.path(dataDirectory, 'data')
  if (!file.exists(data.store)) {
    dir.create(data.store)
  }
  
  browser.store <- file.path(dataDirectory, 'browserFiles')
  if (!file.exists(browser.store)) {
    dir.create(browser.store)
  }
  
  vcf.store <- file.path(dataDirectory, 'VCFfiles')
  if (!file.exists(vcf.store)) {
    dir.create(vcf.store)
  }
  
  ## Loading in list of SNV positions and locations of WC regions
  snvs <- read.table(positions, header=F)
  snvs <- GRanges(seqnames=snvs$V1, IRanges(start=snvs$V2, end=snvs$V2))
  WC.regions <- read.table(WCregions, header=F, sep = ":")
  WC.regions <- GRanges(seqnames=WC.regions$V1, IRanges(start=WC.regions$V2, end=WC.regions$V3), filename=as.character(WC.regions$V4))
  
  ## Run phasing pipeline for selected chromosomes
  for (chr in chromosomes) {
    message("Working on ",chr)
    chr <- as.character(chr) #always consider chromosome name as character
    
    #Select chromosome of interest from the list
    snvs.chr <- snvs[seqnames(snvs) == chr]
    WCregions.chr <- WC.regions[seqnames(WC.regions) == chr]
    seqlevels(snvs.chr) <- chr
    seqlevels(WCregions.chr) <- chr
    
    #load data into matrix  
    matrices <- loadMatrices(bamfilespath=bamfilespath, positions=snvs.chr, WCregions=WCregions.chr, pairedEndReads=pairedEndReads, min.mapq=min.mapq, min.baseq=min.baseq)
    
    #phase data
    srt.matrices <- sortMatrices(matrices, num.iterations=num.iterations)
  
    #filter unreliable data
    assem.haps <- assembleHaps(srt.matrices, translateBases=translateBases, score2qual=score2qual)
    
    #fill gaps in haplotypes
    if (!is.null(fillMissAllele)) {
      assem.haps <- fillGaps(data.object=assem.haps, merged.bam=fillMissAllele, min.mapq=min.mapq, min.baseq=min.baseq, translateBases=translateBases, score2qual=score2qual, chromosome=chr)  
    }
    
    #export phased haplotypes
    destination <- file.path(data.store, paste0(chr, '_phased.RData')) #save original data object with sorted matrices
    save(srt.matrices, file=destination)
    
    destination <- file.path(phased.store, paste0(chr, '_phased_hap1.txt')) #save phased alleles per haplotype in separate files
    write.table(assem.haps$hap1.cons, file=destination, row.names = F)
    destination <- file.path(phased.store, paste0(chr, '_phased_hap2.txt'))
    write.table(assem.haps$hap2.cons, file=destination, row.names = F)
    
    destination <- file.path(phased.store, paste0(chr, '_phasedFiles_hap1.txt'))
    hap1.files <- data.frame(names(assem.haps$hap1.files), do.call(rbind, lapply(assem.haps$hap1.files, rbind)))
    names(hap1.files) <- c("Filenames", "Simil", "Disimil")
    write.table(hap1.files, file=destination, row.names = F)
    destination <- file.path(phased.store, paste0(chr, '_phasedFiles_hap2.txt'))
    hap2.files <- data.frame(names(assem.haps$hap2.files), do.call(rbind, lapply(assem.haps$hap2.files, rbind)))
    names(hap2.files) <- c("Filenames", "Simil", "Disimil")
    write.table(hap2.files, file=destination, row.names = F)
    
    if (!is.null(exportVCF) & !is.null(bsGenome)) {		
    	exportVCF(index = exportVCF, outputDirectory = vcf.store, phasedHap = assem.haps, bsGenome=bsGenome, chromosome = chr)
    }	
  
    #split reads per haplotype  
    if (splitPhasedReads) {
      haps.gr <- splitReads(data.object=assem.haps, bamfilespath=bamfilespath, pairedEndReads=pairedEndReads, min.mapq=min.mapq)
      destination <- file.path(data.store, paste0(chr, '_reads.RData'))
      save(haps.gr, file=destination)
      
      exportBedGraph(index=paste0(chr, '_hap1'), outputDirectory=browser.store, fragments=haps.gr$hap1, col="0,128,255")
      exportBedGraph(index=paste0(chr, '_hap2'), outputDirectory=browser.store, fragments=haps.gr$hap2, col="0,255,255")
    }
    
    #call BreakPointR on phased reads object
    if (callBreaks) {
      breakspath <- file.path(dataDirectory, 'BreakPointR')
      
      if (!file.exists(breakspath)) {
        dir.create(breakspath)
      }
      
      hap1 <- haps.gr$hap1
      hap2 <- haps.gr$hap2
      strand(hap1) <- "-"
      strand(hap2) <- "+"
      phased.haps <- sort(append(hap1,hap2), ignore.strand=T)
      breakpoints <- runBreakpointr(input.data = phased.haps, pairedEndReads=pairedEndReads, chromosomes=chr, windowsize=15000, scaleWindowSize=TRUE, pair2frgm=TRUE)
      writeBedFile(index=chr, outputDirectory=breakspath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
    }
  }
}
