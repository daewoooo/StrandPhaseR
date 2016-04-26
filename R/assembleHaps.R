#' This will take sorted matrices and will calculate concordance of each single cell to the consensus haplotypes in order to assemble highly 
#' accurate haplotypes
#' 
#' @param data.object containing sorted watson and crick haplotypes of each single cell
#' 
#' @author David Porubsky
#' @export

assembleHaps <- function(data.object, translateBases=FALSE, score2qual=FALSE) {
  message("Assembling haplotypes")
  ptm <- proc.time()
  
  hap1 <- data.object[['hap1.bases']]
  hap2 <- data.object[['hap2.bases']]
  hap1.quals <- data.object[['hap1.quals']]
  hap2.quals <- data.object[['hap2.quals']]
  hap1.files <- data.object[['hap1.files']]
  hap2.files <- data.object[['hap2.files']]
  genomic.pos <- data.object[['genomic.pos']]
  
  #Calculate initial consensus haplotypes with default settings
  hap1.cons <- exportConsensus(hap1, hap1.quals)
  hap2.cons <- exportConsensus(hap2, hap2.quals)
  #hap1.cons <- exportConsensus(hap1, hap1.quals, filt.ambig=T, min.cov=1)
  #hap2.cons <- exportConsensus(hap2, hap2.quals, filt.ambig=T, min.cov=1)
  
  ## obtain HET positions only
  cov.both <- intersect(hap1.cons$pos, hap2.cons$pos)
  hap1.both <- hap1.cons[hap1.cons$pos %in% cov.both,]
  hap2.both <- hap2.cons[hap2.cons$pos %in% cov.both,]
  HETpos <- hap1.both[hap1.both$bases != hap2.both$bases,]$pos
  
  assembled.hap1 <- list()
  assembled.hap2 <- list()
  mask.hap1 <- vector()
  mask.hap2 <- vector()
  for (i in 1:length(hap1.files)) {
    #message("Iteration ",i)
    
    ## cell hap1
    hap1.cell <- data.frame(pos=1:length(hap1[i,]), bases=hap1[i,])
    hap1.cell <- hap1.cell[hap1.cell$base > 0,]
    hap1.cell.HET <- hap1.cell[hap1.cell$pos %in% HETpos,]
    hap1.cell.bases <- hap1.cell.HET$bases
    
    ## Consensus bases
    hap1.consBases <- hap1.cons[hap1.cons$pos %in% hap1.cell.HET$pos,]$bases
    hap2.consBases <- hap2.cons[hap2.cons$pos %in% hap1.cell.HET$pos,]$bases
    
    ## Compare hap1 cell to consensus hap1 and hap1 and calculate concordance as percentage of matched allele in comparison to consensus
    match.cis <- length(hap1.cell.bases[hap1.cell.bases == hap1.consBases])
    match.trans <- length(hap1.cell.bases[hap1.cell.bases == hap2.consBases])
    hap1.cis.concordance <- match.cis / length(hap1.cell.bases)
    hap1.trans.concordance <- match.trans / length(hap1.cell.bases)
    
    ## cell hap2
    hap2.cell <- data.frame(pos=1:length(hap2[i,]), bases=hap2[i,])
    hap2.cell <- hap2.cell[hap2.cell$base > 0,]
    hap2.cell.HET <- hap2.cell[hap2.cell$pos %in% HETpos,]
    hap2.cell.bases <- hap2.cell.HET$bases
    
    ## Consensus bases
    hap1.consBases <- hap1.cons[hap1.cons$pos %in% hap2.cell.HET$pos,]$bases
    hap2.consBases <- hap2.cons[hap2.cons$pos %in% hap2.cell.HET$pos,]$bases
    
    ## Compare hap2 cell to consensus hap1 and hap1 and calculate concordance as percentage of matched allele in comparison to consensus
    match.cis <- length(hap2.cell.bases[hap2.cell.bases == hap2.consBases])
    match.trans <- length(hap2.cell.bases[hap2.cell.bases == hap1.consBases])
    hap2.cis.concordance <- match.cis / length(hap2.cell.bases)
    hap2.trans.concordance <- match.trans / length(hap2.cell.bases)
    
    ## Filter out unreliably phased cells
    diff.level.hap1 <- ((abs(hap1.cis.concordance - hap1.trans.concordance))/(hap1.cis.concordance + hap1.trans.concordance)/2)*100
    diff.level.hap2 <- ((abs(hap2.cis.concordance - hap2.trans.concordance))/(hap2.cis.concordance + hap2.trans.concordance)/2)*100
    
    if (diff.level.hap1>25 & hap1.cis.concordance>0.9 & nrow(hap1.cell.HET)>0) {
      hap1.file <- hap1.files[i]
      assembled.hap1[[hap1.file]] <- c(hap1.cis.concordance, hap1.trans.concordance)  
    } else {
      mask.hap1 <- c(mask.hap1, i)
    }
    
    if (diff.level.hap2>25 & hap2.cis.concordance>0.9 & nrow(hap2.cell.HET)>0) {
      hap2.file <- hap2.files[i]
      assembled.hap2[[hap2.file]] <- c(hap2.cis.concordance, hap2.trans.concordance)  
    } else {
      mask.hap2 <- c(mask.hap2, i)
    }
  }

  ## Remove unreliable single cell haplotypes
  if (length(mask.hap1)>0) {
    hap1.files <- hap1.files[-mask.hap1]
    hap1 <- hap1[-mask.hap1,]
    hap1.quals <- hap1.quals[-mask.hap1,]
  }  
  
  if (length(mask.hap2)>0) {
    hap2.files <- hap2.files[-mask.hap2]
    hap2 <- hap2[-mask.hap2,]
    hap2.quals <- hap2.quals[-mask.hap2,]
  }
  
  ## Recalculate consensus haplotyes after filtering unreliable single cell haplotypes (using user defined settings)
  hap1.cons <- exportConsensus(hap1, hap1.quals, filt.ambig=FALSE, min.cov=1, translateBases=translateBases, score2qual=score2qual)
  hap2.cons <- exportConsensus(hap2, hap2.quals, filt.ambig=FALSE, min.cov=1, translateBases=translateBases, score2qual=score2qual)
  
  hap1.GenomicPos <- genomic.pos[hap1.cons$pos]
  hap2.GenomicPos <- genomic.pos[hap2.cons$pos]
  
  hap1.cons$pos <- hap1.GenomicPos
  hap2.cons$pos <- hap2.GenomicPos
  
  ## Export data as a single data object (List)
  assem.haps <- list()
  #assem.haps[['hap1.bases']] <- hap1
  #assem.haps[['hap2.bases']] <- hap2
  #assem.haps[['hap1.quals']] <- hap1.quals
  #assem.haps[['hap2.quals']] <- hap2.quals
  assem.haps[['hap1.cons']] <- hap1.cons
  assem.haps[['hap2.cons']] <- hap2.cons
  assem.haps[['hap1.files']] <- assembled.hap1
  assem.haps[['hap2.files']] <- assembled.hap2
  
  time <- proc.time() - ptm
  message("Time spent: ",round(time[3],2),"s")
  
  return(assem.haps)
}    