#' This function uses circular-binary-segementation to segment single cell haplotypes based on their similarity to consensus haploypes of all cells.
#' 
#' @param data.object A \code{list} output from the function \code{\link{compareSingleCellHaps}}
#' @param chromosome  A chromosome id.
#' @importFrom fastseg fastseg
#' @importFrom IRanges IRanges
#' @author David Porubsky
#' @export

LOHseeker <- function(data.object=NULL, chromosome=NULL, bin.size=10) {
  message(" Searching for LOH", appendLF=F); ptm <- proc.time()
  
  all.segments <- GRangesList()
  for (i in 1:length(data.object)) {
    #get data
    cell.data <- data.object[[i]]
    cell.haps <- split(cell.data, cell.data$CellID)
    if (length(cell.haps) < 2) {next}
    cell.comparison.H1 <- cell.haps[[1]]
    cell.comparison.H2 <- cell.haps[[2]]
    if (nrow(cell.comparison.H1) < bin.size | nrow(cell.comparison.H2) < bin.size) {next}
    #get cell ID
    cell.comparison.H1.id <- names(cell.haps[1])
    cell.id <- gsub(cell.comparison.H1.id, pattern = "_H1", replacement = "")
    
    #segmentation
    H1.segs <- fastseg::fastseg(cell.comparison.H1$cons1.simil, segMedianT = 0.5) #comparing only to one consesus hap because the other is the complete mirror
    H2.segs <- fastseg::fastseg(cell.comparison.H2$cons1.simil, segMedianT = 0.5)
    
    #filter segments with mean similarity values in between 0.1 and 0.9 
    H1.segs <- H1.segs[H1.segs$seg.mean < 0.1 | H1.segs$seg.mean > 0.9]
    H2.segs <- H2.segs[H2.segs$seg.mean < 0.1 | H2.segs$seg.mean > 0.9]
    
    if (length(H1.segs)>0 & length(H2.segs)>0) {
      #get genomic positions of localized segments
      H1.segs.ranges <- ranges(H1.segs)
      H2.segs.ranges <- ranges(H2.segs)
      H1.starts <- cell.comparison.H1$start[start(H1.segs.ranges)]
      H1.ends <- cell.comparison.H1$end[end(H1.segs.ranges)]
      H2.starts <- cell.comparison.H2$start[start(H2.segs.ranges)]
      H2.ends <- cell.comparison.H2$end[end(H2.segs.ranges)]
      
      #get haplotype for each segment (simil close to 1 means that given haplotype belongs to compared consesusHap1 or vice versa)
      H1.segs.haps <- ifelse(H1.segs$seg.mean > 0.9, "H1", ifelse(H1.segs$seg.mean < 0.1, 'H2', 'NA'))
      H2.segs.haps <- ifelse(H2.segs$seg.mean > 0.9, "H1", ifelse(H2.segs$seg.mean < 0.1, 'H2', 'NA'))
      
      #create GRanges object of loacalized segments
      H1.segments <- GRanges(seqnames=chromosome, ranges=IRanges::IRanges(start=H1.starts, end=H1.ends), exp.Hap="H1", obs.Hap=H1.segs.haps, ID=cell.id)
      H2.segments <- GRanges(seqnames=chromosome, ranges=IRanges::IRanges(start=H2.starts, end=H2.ends), exp.Hap="H2", obs.Hap=H2.segs.haps, ID=cell.id)
      cell.segments <- c(H1.segments, H2.segments)
      all.segments[[cell.id]] <- cell.segments
    }  
  }
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  if (length(all.segments) > 0) {
    return(unlist(all.segments))
  } else {
    return(NULL) 
  }
}