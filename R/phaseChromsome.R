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
#' @param inputfolder Path to the bam files to process
#' @param outputfolder Output directory. If non-existent it will be created.
#' @param positions Filename with listed position of SNVs for given chromosome (format: chrName SNVpos).
#' @param WCregions Filename of all WC region for a given chromosome (format: chrName:Start:End:FileName).
#' @param chromosome If only a subset of the chromosomes should be processed, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param min.baseq Minimum base quality to consider a base for phasing.
#' @param num.iterations Number of iteration to sort watson and crick matrices.
#' @param translateBases ...
#' @param fillGaps ...
#' @param splitPhasedReads Set to \code{TRUE} if you want to split reads per haplotype.
#' @param compareSingleCells Set to \code{TRUE} if you want to compare haplotypes at the single-cell level.
#' @param callBreaks ...
#' @param exportVCF ...
#' @param bsGenome A \code{BSgenome} object which contains the refernce DNA sequence

#' @author David Porubsky
#' @export

phaseChromosome <- function(inputfolder, outputfolder='./StrandPhaseR_analysis', positions=NULL, WCregions=NULL, chromosome=NULL, pairedEndReads=TRUE, min.mapq=10, min.baseq=30, num.iterations=2, translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, callBreaks=FALSE, exportVCF=NULL, bsGenome=NULL) {
  
  message("Working on chromosome ",chromosome)
  
  ## Check user input
  if (!is.null(exportVCF) & is.null(bsGenome)) {
    warning("VCf file cannot be created because reference genome is NULL")
  }	
  
  phased.store <- file.path(outputfolder, 'Phased')
  data.store <- file.path(outputfolder, 'data')
  browser.store <- file.path(outputfolder, 'browserFiles')
  vcf.store <- file.path(outputfolder, 'VCFfiles')
  singlecell.store <- file.path(outputfolder, 'SingleCellHaps')
  
  #load data into matrix
  matrices <- loadMatrices(inputfolder=inputfolder, positions=positions, WCregions=WCregions, pairedEndReads=pairedEndReads, min.mapq=min.mapq, min.baseq=min.baseq)
  
  #Check if sufficient data were loaded
  if (length(matrices) > 0) {
  
    #phase data
    srt.matrices <- sortMatrices(matrices, num.iterations=num.iterations)
  
    #filter unreliable data
    assem.haps <- assembleHaps(srt.matrices, translateBases=translateBases)
    
    #fill gaps in haplotypes
    header <- read.table(fillMissAllele, stringsAsFactors = FALSE, fill=TRUE, comment.char = "&", nrows = 1)
    if (grepl(header, pattern = "VCF", ignore.case = TRUE)) {
      assem.haps <- fillGapsWithVCF(data.object=assem.haps, ref.vcf=fillMissAllele, chromosome=chromosome)
    }
    
    if (grepl(fillMissAllele, pattern = "\\.bam$")) {
      assem.haps <- fillGapsWithBam(data.object=assem.haps, merged.bam=fillMissAllele, min.mapq=min.mapq, min.baseq=min.baseq, translateBases=translateBases, chromosome=chromosome)  
    }
  
    #compara single-cell haplotypes to assembled consensus haplotypes
    if (compareSingleCells) {
      #compare single cell haplotypes to the consensus haplotypes
      suppressWarnings( cell.comparisons.l <- compareSingleCellHaps(consensusHaps=assem.haps, sortedHaps=srt.matrices) )
      #plt.df <- melt(cell.comparisons.l, measure.vars = c('cons1.simil','cons2.simil'))  
      #plt <- ggplot(plt.df, aes(y=value,x=start, color=variable)) + geom_step()  + facet_grid(CellID ~ .) + theme_bw() + theme(strip.text.y = element_text(angle=0), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + scale_color_manual(values = c("darkgoldenrod1", "dodgerblue2"))
      #destination <- file.path(singlecell.store, paste0(chromosome, '_singleCellHaps.pdf'))
      #suppressMessages( ggsave(destination, plot=plt, device="pdf", width=10, height=length(unique(plt.df$CellID))*0.5, limitsize=F) )
  
      if (!is.null(cell.comparisons.l)) {
        destination <- file.path(singlecell.store, paste0(chromosome, '_singleCellHaps.pdf'))
        suppressWarnings( plotSingleCellHaps(data=cell.comparisons.l, file=destination) )
        
        #Detect LOH regions
        LOH.regions <- LOHseeker(data.object=cell.comparisons.l, chromosome=chromosome)
        LOH.regions.df <- data.frame(LOH.regions)
        destination <- file.path(singlecell.store, paste0(chromosome, '_singleCell_LOH.txt'))
        write.table(LOH.regions.df, file = destination, quote = F, row.names = F)
      }  
    }  
  
    #add chromosome name to exported phased files 
    chrName.hap1 <- data.frame(chr=rep(chromosome, nrow(assem.haps$hap1.cons)))
    assem.haps$hap1.cons <- cbind(chrName.hap1, assem.haps$hap1.cons)
    chrName.hap2 <- data.frame(chr=rep(chromosome, nrow(assem.haps$hap2.cons)))
    assem.haps$hap2.cons <- cbind(chrName.hap2, assem.haps$hap2.cons)	 	
  
    #export phased haplotypes
    destination <- file.path(data.store, paste0(chromosome, '_phased.RData')) #save original data object with sorted matrices
    save(srt.matrices, file=destination)
    
    destination <- file.path(phased.store, paste0(chromosome, '_phased_hap1.txt')) #save phased alleles per haplotype in separate files
    write.table(assem.haps$hap1.cons, file=destination, row.names = F)
    destination <- file.path(phased.store, paste0(chromosome, '_phased_hap2.txt'))
    write.table(assem.haps$hap2.cons, file=destination, row.names = F)
    
    destination <- file.path(phased.store, paste0(chromosome, '_phasedFiles_hap1.txt'))
    hap1.files <- data.frame(names(assem.haps$hap1.files), do.call(rbind, lapply(assem.haps$hap1.files, rbind)))
    names(hap1.files) <- c("Filenames", "Simil", "Disimil")
    write.table(hap1.files, file=destination, row.names = F)
    destination <- file.path(phased.store, paste0(chromosome, '_phasedFiles_hap2.txt'))
    hap2.files <- data.frame(names(assem.haps$hap2.files), do.call(rbind, lapply(assem.haps$hap2.files, rbind)))
    names(hap2.files) <- c("Filenames", "Simil", "Disimil")
    write.table(hap2.files, file=destination, row.names = F)
    destination <- file.path(phased.store, 'phased_haps.txt')
    write.table(data.frame(assem.haps$assem.haps), file=destination, row.names = F, col.names=F, quote = F, append = T, sep="\t")	  

    if (!is.null(exportVCF) & !is.null(bsGenome)) {		
      exportVCF(index = exportVCF, outputfolder = vcf.store, phasedHap = assem.haps, bsGenome=bsGenome, chromosome=chromosome)
    }	
    
    #split reads per haplotype  
    if (splitPhasedReads) {
      haps.gr <- splitReads(data.object=assem.haps, inputfolder=inputfolder, pairedEndReads=pairedEndReads, min.mapq=0, filterAltAlign=FALSE)
      destination <- file.path(data.store, paste0(chromosome, '_reads.RData'))
      save(haps.gr, file=destination)
      
      exportBedGraph(index=paste0(chromosome, '_hap1'), outputfolder=browser.store, fragments=haps.gr$hap1, col="0,128,255")
      exportBedGraph(index=paste0(chromosome, '_hap2'), outputfolder=browser.store, fragments=haps.gr$hap2, col="0,255,255")
    }
  
  } else {
    message(" Insufficient data to assemble haplotypes, skipping ...")
  }  
  
  #call BreakPointR on phased reads object
#  if (callBreaks) {
#    breakspath <- file.path(outputfolder, 'BreakPointR')
      
#    if (!file.exists(breakspath)) {
#      dir.create(breakspath)
#    }
      
#    hap1 <- haps.gr$hap1
#    hap2 <- haps.gr$hap2
#    strand(hap1) <- "-"
#    strand(hap2) <- "+"
#    phased.haps <- sort(append(hap1,hap2), ignore.strand=T)
#   breakpoints <- runBreakpointr(bamfile = phased.haps, pairedEndReads=pairedEndReads, chromosomes=chromosome, windowsize=50, binMethod="reads", pair2frgm=TRUE)	
#    writeBedFile(index=chromosome, outputDirectory=breakspath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
#  }

} #end of function
