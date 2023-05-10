#' Map meiotic recombination events in a single family trio [father-mother-child].
#' 
#' This function takes as an input phased VCF files for each member of a family trio and maps breakpoints
#' of meiotic recombination in each inherited parental homolog.
#' 
#' @param parent1 A path to a VCF file to be loaded for a parent 1.
#' @param parent2 A path to a VCF file to be loaded for a parent 2.
#' @param child A path to a VCF file to be loaded for a child.
#' @param method A user defined method to be used to map changes in haplotype blocks [default: CBS]
#' @inheritParams findRecomb
#' @inheritParams vcf2vranges
#' @importFrom fastseg fastseg
#' @importFrom dplyr recode
#' @return A \code{list} object that contains mapped meiotic breakpoints and inherited hapltype segments 
#' for each homolog in a child.
#' @author David Porubsky
#' @export
#' 
mapRecomb <- function(parent1=NULL, parent2=NULL, child=NULL, genome='hg38', method='CBS', minSeg=100, smooth=3, collapse.amb=TRUE) {
  
  ## Load VCF files
  if (!is.null(parent1) & file.exists(parent1)) {
    par1 <- vcf2vranges(vcfFile = parent1, genoField = 'GT', translateBases = TRUE, genome = genome)
  } else {
    ## Create a dummy VCF record
    par1 <- VRanges(seqnames = 'chr', ranges = IRanges(start=1, end=1), 
                    ref = 'A', alt = 'N', 
                    totalDepth = 0, refDepth = 0, altDepth = 0, 
                    sampleNames = 'parent1', softFilterMatrix = matrix(FALSE))
  }
  if (!is.null(parent2) & file.exists(parent2)) {
    par2 <- vcf2vranges(vcfFile = parent2, genoField = 'GT', translateBases = TRUE, genome = genome)
  } else {
    ## Create a dummy VCF record
    par2 <- VRanges(seqnames = 'chr', ranges = IRanges(start=1, end=1), 
                    ref = 'A', alt = 'N', 
                    totalDepth = 0, refDepth = 0, altDepth = 0, 
                    sampleNames = 'parent2', softFilterMatrix = matrix(FALSE))
  }
  if (!is.null(child) & file.exists(child)) {
    child <- vcf2vranges(vcfFile = child, genoField = 'GT', translateBases = TRUE, genome = genome)
  }
  
  ## Keep only child's HET SNVs
  mask <- child$H1 != 'N' & child$H2 != 'N' & child$H1 == child$H2
  child <- child[!mask]
  
  ## Get only shared SNVs between datasets
  comparison.obj <- child
  comparison.obj$par1.H1 <- 'N'
  comparison.obj$par1.H2 <- 'N'
  comparison.obj$par2.H1 <- 'N'
  comparison.obj$par2.H2 <- 'N'
  if (!is.null(parent1) & file.exists(parent1)) {
    shared.par1 <- findOverlaps(par1, child)
    comparison.obj$par1.H1[subjectHits(shared.par1)] <- par1$H1[queryHits(shared.par1)]
    comparison.obj$par1.H2[subjectHits(shared.par1)] <- par1$H2[queryHits(shared.par1)]
  }
  if (!is.null(parent2) & file.exists(parent2)) {
    shared.par2 <- findOverlaps(par2, child)
    comparison.obj$par2.H1[subjectHits(shared.par2)] <- par2$H1[queryHits(shared.par2)]
    comparison.obj$par2.H2[subjectHits(shared.par2)] <- par2$H2[queryHits(shared.par2)]
  }  
  
  ## Compare child.H1 to par1
  comparison.obj$c1_to_par1 <- 'N'
  mask <- comparison.obj$H1 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2
  comparison.obj$c1_to_par1[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par1.H1[mask], 'P1', 'P2')
  c1_to_par1.breaks <- length(rle(comparison.obj$c1_to_par1[mask])$lengths)
  ## Compare child.H2 to par1
  comparison.obj$c2_to_par1 <- 'N'
  mask <- comparison.obj$H2 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2
  comparison.obj$c2_to_par1[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par1.H1[mask], 'P1', 'P2')
  c2_to_par1.breaks <- length(rle(comparison.obj$c2_to_par1[mask])$lengths)
  ## Compare child.H1 to par2
  comparison.obj$c1_to_par2 <- 'N'
  mask <- comparison.obj$H1 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2
  comparison.obj$c1_to_par2[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par2.H1[mask], 'P1', 'P2')
  c1_to_par2.breaks <- length(rle(comparison.obj$c1_to_par2[mask])$lengths)
  ## Compare child.H2 to par2
  comparison.obj$c2_to_par2 <- 'N'
  mask <- comparison.obj$H2 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2
  comparison.obj$c2_to_par2[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par2.H1[mask], 'P1', 'P2')
  c2_to_par2.breaks <- length(rle(comparison.obj$c2_to_par2[mask])$lengths)
  ## Get parental homologs
  if ((c1_to_par1.breaks + c2_to_par2.breaks) < (c2_to_par1.breaks + c1_to_par2.breaks)) {
    comparisons <- comparison.obj[,c('c1_to_par1', 'c2_to_par2')]
    inherited.homologs <- c(paste0('H1.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par1))),
                            paste0('H2.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par2)))
    )
    H1.inherited <- as.character(runValue(sampleNames(par1)))
    H2.inherited <- as.character(runValue(sampleNames(par2)))
  } else {
    comparisons <- comparison.obj[,c('c1_to_par2', 'c2_to_par1')]
    inherited.homologs <- c(paste0('H1.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par2))),
                            paste0('H2.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par1)))
    )
    H1.inherited <- as.character(runValue(sampleNames(par2)))
    H2.inherited <- as.character(runValue(sampleNames(par1)))
  }
  names(mcols(comparisons)) <- c('H1.comp', 'H2.comp')
  
  if (method == 'CBS') {
    H1.recomb <- findRecomb(comparison = comparisons[,'H1.comp'], minSeg = minSeg, smooth = smooth, collapse.amb = collapse.amb)
    H2.recomb <- findRecomb(comparison = comparisons[,'H2.comp'], minSeg = minSeg, smooth = smooth, collapse.amb = collapse.amb)
  }
  
  ## Add field of inherited homologs
  H1.recomb$hap.segm$inherited <- H1.inherited
  H1.recomb$recomb.break$inherited <- H1.inherited
  H2.recomb$hap.segm$inherited <- H2.inherited
  H2.recomb$recomb.break$inherited <- H2.inherited
  
  ## Return results object
  return(list(comparisons = comparisons, 
              H1.segments = H1.recomb$hap.segm, 
              H1.recomb.breaks = H1.recomb$recomb.break, 
              H2.segments = H2.recomb$hap.segm,
              H2.recomb.breaks = H2.recomb$recomb.break,
              inherited.homologs = inherited.homologs)
  )
} 


#' Find meiotic recombination breakpoints using circular binary segmentation (CBS)
#' 
#' This function takes as an input a \code{\link{VRanges-class}} object with extra metacolumn that contains
#' for each variant parantal identity. P1 -> parent1 & P2 -> parent2.
#' 
#' @param comparison A \code{\link{VRanges-class}} object 
#' @param minSeg Minimal length (number of variants) being reported as haplotype block (\code{fastseg} parameter).
#' @param smooth Number of consecutive variants being considered as a random error and so being corrected (flipped).
#' @param collapse.amb Set to \code{TRUE} if segments with ambiguous haplotype assignments should be collapsed.
#' @importFrom fastseg fastseg
#' @importFrom dplyr recode
#' @author David Porubsky
#' @export

findRecomb <- function(comparison=NULL, minSeg=100, smooth=3, collapse.amb=TRUE) {
  
  ## Helper function
  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1  
    } 
  }
  
  collapseBins <- function(gr, id.field=0) {
    ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
  }
  
  ## Remove missing values
  comparison.filt <- comparison[mcols(comparison)[,1] != 'N']
  
  if (length(comparison.filt) >= 2*minSeg) {
    ## Recode comparison into in 0/1 vector
    comparison.filt$comp.vector <- dplyr::recode(mcols(comparison.filt)[,1], 'P1' = 0, 'P2' = 1)
    
    ## Run CBS segmentation on 0/1 vector 
    segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    
    while (any(segs$num.mark <= smooth)) {
      toSwitch <- which(segs$num.mark <= smooth)
      switch.segs <- segs[toSwitch]
      switch.pos <- mapply(function(x,y) {x:y}, x=switch.segs$startRow, y=switch.segs$endRow)
      switch.pos <- unlist(switch.pos)
      
      switched.vals <- sapply(comparison.filt$comp.vector[switch.pos], switchValue) #SWITCH
      #comparison$comp.vector <- comparison$comp.vector[-switch.pos]  #DELETE
      comparison.filt$comp.vector[switch.pos] <- switched.vals
      segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    }
    gen.ranges <- IRanges(start=start(comparison.filt)[segs$startRow], end=end(comparison.filt)[segs$endRow])
    ranges(segs) <- gen.ranges
    ## Add chromosome name
    seqlevels(segs) <- seqlevels(comparison.filt)
    #segs$index <- index
    
    segs$match[segs$seg.mean <= 0.25] <- 'hap1'
    segs$match[segs$seg.mean >= 0.75] <- 'hap2'
    segs$match[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'amb'
    
    ## Remove segments with mixed H1 and H2 signal
    if (collapse.amb) {
      segs <- segs[segs$match != 'amb']
    }
    
  } else {
    message("    Low density of informative SNVs, skipping ...")
    segs <- GenomicRanges::GRanges(seqnames=seqlevels(comparison), 
                                   ranges=IRanges(start=min(start(comparison)), end=max(end(comparison))),
                                   ID=NA, num.mark=0, seg.mean=0, startRow=0, endRow=0, match=NA
    )
  }
  
  ## Get haplotype segments
  if (length(segs) > 1) {
    segm <- collapseBins(gr = segs, id.field = 6)
  } else {
    segm <- segs
  }
  
  ## Meiotic recombination breakpoints
  if (length(segm) > 1) {
    suppressWarnings( recomb.break <- gaps(segm, start = start(segm)) )
    start(recomb.break) <- start(recomb.break) - 1
    end(recomb.break) <- end(recomb.break) + 1
  } else {
    ## Create a dummy breakpoint
    recomb.break <- GenomicRanges::GRanges(seqnames = seqlevels(comparison), ranges = IRanges(start=1, end=1))
  } 
  return(list(hap.segm = segm, recomb.break = recomb.break))
}