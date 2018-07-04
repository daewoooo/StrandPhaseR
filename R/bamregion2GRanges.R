#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param region If only a subset of the genomic regions should be loaded.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag testPairedEndBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author Aaron Taudt, David Porubsky, Ashley Sanders
#' @export

bamregion2GRanges <- function(bamfile, bamindex=bamfile, region=NULL, pairedEndReads=FALSE, min.mapq=10, filterAltAlign=TRUE) {

  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  
  ## Check if bam is truly paired ended in case pairedEndReads set to TRUE
  is.Paired <- Rsamtools::testPairedEndBam(file = bamfile, index = bamindex)
  if (pairedEndReads) {
    if (!is.Paired) {
      warning("You are trying to process single-ended BAM as paired-ended, Please set proper BAM directioanlity!!!")
    } 
  } else {
    if (is.Paired) {
      warning("You are trying to process paired-ended BAM as single-ended, Please set proper BAM directioanlity!!!")
    }  
  }
  
  ## read in reads data
  if (pairedEndReads) {
    #suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(region), what=c('seq', 'qual','mapq','cigar'), flag=scanBamFlag(isDuplicate=F))) )
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(region), what=c('seq', 'qual','mapq','cigar','flag'), )) )	
  } else {
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(region), what=c('seq', 'qual','mapq','cigar'), flag=scanBamFlag(isDuplicate=F))) )
  } 
  
  ## Second mate of the pair will inherit directionality from the first mate of the pair
  if (pairedEndReads) {
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
    strand(data.last) <- strand(data.first)
    data <- GenomicRanges::sort(c(data.first, data.last), ignore.strand=TRUE)
  } else {
    data <- as(data.raw, 'GRanges')
  }
  
  ## Filter duplicates for pairedEndReads
  if (pairedEndReads) {
    bit.flag <- bitwAnd(1024, data$flag)
    mask <- bit.flag == 0 	
    data <- data[mask]
  }  

  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }
  
  ## filter XA tag
  if (filterAltAlign) {
    data <- data[is.na(mcols(data)$XA)]
  }    
  
  #data <- data[seqnames(data) %in% seqlevels(region)]
  #seqlevels(data) <- seqlevels(region)
  data <- GenomeInfoDb::keepSeqlevels(data, seqlevels(region), pruning.mode="coarse")
  #data <- GenomeInfoDb::keepSeqlevels(data, seqlevels(region))	
  	
  return(data)
}
