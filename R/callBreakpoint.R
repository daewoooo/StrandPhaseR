#' This function will call BreakPointR on phased reads object
#' 
#' @param datapath
#' @param outputfolder
#' @import breakpointR
#' 
#' @author David Porubsky

callBreakpointR <- function(datapath, outputfolder='BreakPointR') {
  
  breakspath <- file.path(datapath, outputfolder)
  
  if (!file.exists(breakspath)) {
    dir.create(breakspath)
  }
  
  RData.objects <- list.files(datapath, pattern="_reads.RData$", full=T)
  
  data <- get(load(RData.objects)) 
  
  hap1 <- data$hap1
  hap2 <- data$hap2
  strand(hap1) <- "-"
  strand(hap2) <- "+"
  #hap1.dis <- disjoin(hap1)
  #hap2.dis <- disjoin(hap2)
  
  #dup.hap1 <- mapply(function(X,Y) all(X,Y), X=duplicated(start(hap1)), Y=duplicated(end(hap1)))
  #dup.hap2 <- mapply(function(X,Y) all(X,Y), X=duplicated(start(hap2)), Y=duplicated(end(hap2)))
  #hap1 <- hap1[!dup.hap1]
  #hap2 <- hap2[!dup.hap2]

  phased.haps <- sort(append(hap1,hap2), ignore.strand=T)
  
  #phased.haps <- phased.haps[width(phased.haps) > mean(width(phased.haps))]
  #phased.haps <- phased.haps[width(phased.haps) >= 76]
  #phased.haps$mapq <- 40
  
  breakpoints <- runBreakpointr(input.data = phased.haps, pairedEndReads = F, chromosomes = "chr22", windowsize = 15000, scaleWindowSize = T)
  writeBedFile(index=outputfolder, outputDirectory=breakspath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
}  