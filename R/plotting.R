#' This function will take phased reads for each chromosome and plot them 
#' 
#' @param data.object containing split watson and crick reads per haplotype
#' @param plotChromosomes
#' @param file name of the file to save produced plots
#' 
#' @author David Porubsky
#'

plotPhasedReads <- function(datapath, plotChromosomes = NULL, file=NULL) {
  
  RData.objects <- list.files(datapath, pattern="_reads.RData$", full=T)
  
  data <- get(load(RData.objects))  
  
  hap1 <- data$hap1
  hap2 <- data$hap2
  
  #get coverage per genomic position
  rle.hap1.cov <- coverage(hap1)
  hap1.cov <- runValue(rle.hap1.cov) # get values of coverage per genomic site 
  hap1.cov <- unlist(hap1.cov, use.names = FALSE) # transfer coverage values into vector
  hap1.cov.ranges <- ranges(rle.hap1.cov) # get genomic ranges with certain value of coverage 
  hap1.cov.ranges <- unlist(hap1.cov.ranges, use.names = FALSE) 
  df.hap1Cov <- data.frame(start=start(hap1.cov.ranges), end=end(hap1.cov.ranges), cov=hap1.cov)
  
  rle.hap2.cov <- coverage(hap2)
  hap2.cov <- runValue(rle.hap2.cov)
  hap2.cov <- unlist(hap2.cov, use.names = FALSE)
  hap2.cov.ranges <- ranges(rle.hap2.cov)
  hap2.cov.ranges <- unlist(hap2.cov.ranges, use.names = FALSE)
  df.hap2Cov <- data.frame(start=start(hap2.cov.ranges), end=end(hap2.cov.ranges), cov=hap2.cov)
  
  ggplt <- ggplot(df.hap1Cov, aes(x=start, y=cov))
  ggplt <- ggplt + geom_step(col='red')
  ggplt <- ggplt + geom_step(data=df.hap2Cov, aes(x=start, y=-cov), col = 'blue', inherit.aes = F)
  
  return(ggplt)
}
    
