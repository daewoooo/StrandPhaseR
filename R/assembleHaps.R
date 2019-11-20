#' This will take sorted matrices and will calculate concordance of each single cell to the consensus haplotypes in order to assemble highly 
#' accurate haplotypes
#' 
#' @param data.object containing sorted watson and crick haplotypes of each single cell
#' @param concordance Level of agreement between single cell and consensus haplotypes
#' @inheritParams exportConsensus
#' 
#' @author David Porubsky
#' @export

assembleHaps <- function(data.object, translateBases=FALSE, concordance=0.9) {
  message(" Assembling haplotypes", appendLF=F); ptm <- proc.time()
  
  hap1 <- data.object[['hap1.bases']]
  hap2 <- data.object[['hap2.bases']]
  hap1.quals <- data.object[['hap1.quals']]
  hap2.quals <- data.object[['hap2.quals']]
  hap1.files <- data.object[['hap1.files']]
  hap2.files <- data.object[['hap2.files']]
  genomic.pos <- data.object[['genomic.pos']]
  
  #Calculate initial consensus haplotypes with default settings
  hap1.cons <- exportConsensus(hap1, hap1.quals, min.cov=2)
  hap2.cons <- exportConsensus(hap2, hap2.quals, min.cov=2)
  
  #if any consensus sequence was return as zero recall consensus with min.cov=1
  if (is.null(nrow(hap1.cons)) | is.null(nrow(hap2.cons))) {
    hap1.cons <- exportConsensus(hap1, hap1.quals, min.cov=1)
    hap2.cons <- exportConsensus(hap2, hap2.quals, min.cov=1)
  }
  consensus.density <- mean(c(nrow(hap1.cons), nrow(hap2.cons)))
  
  ## obtain HET positions only
  cov.both <- intersect(hap1.cons$pos, hap2.cons$pos)
  hap1.both <- hap1.cons[hap1.cons$pos %in% cov.both,]
  hap2.both <- hap2.cons[hap2.cons$pos %in% cov.both,]
  HETpos <- hap1.both[hap1.both$bases != hap2.both$bases,]$pos
  
  #check if the number heterozygous snvs is at least 10% of consensus SNVs from both haplotypes
  #and there is at least 5 HET positions in total
  if (length(HETpos) > 0.1*consensus.density & length(HETpos) >= 5) {
  #if (length(HETpos) >= 10) { #at least 10 overlapping HET positions  
  
    assembled.hap1 <- list()
    assembled.hap2 <- list()
    assembled.haps <- list()	
    mask.hap1 <- vector()
    mask.hap2 <- vector()
    for (i in 1:length(hap1.files)) {
    
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
    
      ## export all phased regions (unfiltered)
      phased.info <- unlist(strsplit(hap1.files[i], '__'))
      ID <- phased.info[1]	
      hap1.phase <- phased.info[3]	
      phased.region <-  phased.info[2]
      phased.region <- unlist(strsplit(phased.region, ':|-'))		
      hap2.phase <- unlist(strsplit(hap2.files[i], '__'))[3]
      phase <- paste0(hap1.phase, hap2.phase)
      assembled.haps[[i]] <- paste(sep = "\t", c("NA", ID, phased.region, phase, hap1.cis.concordance, hap1.trans.concordance, hap2.cis.concordance, hap2.trans.concordance))		

      ## Filter out unreliably phased cells
      diff.level.hap1 <- ((abs(hap1.cis.concordance - hap1.trans.concordance))/(hap1.cis.concordance + hap1.trans.concordance)/2)*100
      diff.level.hap2 <- ((abs(hap2.cis.concordance - hap2.trans.concordance))/(hap2.cis.concordance + hap2.trans.concordance)/2)*100
    
      hap1.file <- hap1.files[i]
      assembled.hap1[[hap1.file]] <- c(hap1.cis.concordance, hap1.trans.concordance)  
      hap2.file <- hap2.files[i]
      assembled.hap2[[hap2.file]] <- c(hap2.cis.concordance, hap2.trans.concordance)  
      
      if (diff.level.hap1 < 25 | hap1.cis.concordance < concordance | nrow(hap1.cell.HET) == 0) {
      #if (!diff.level.hap1 > 25 & hap1.cis.concordance > concordance & nrow(hap1.cell.HET) > 0) {
      #  hap1.file <- hap1.files[i]
      #  assembled.hap1[[hap1.file]] <- c(hap1.cis.concordance, hap1.trans.concordance)  
      #} else {
        mask.hap1 <- c(mask.hap1, i)
      }
    
      #if (diff.level.hap2 > 25 & hap2.cis.concordance > concordance & nrow(hap2.cell.HET) > 0) {
      if (diff.level.hap2 < 25 | hap2.cis.concordance < concordance | nrow(hap2.cell.HET) == 0) {
      #  hap2.file <- hap2.files[i]
      #  assembled.hap2[[hap2.file]] <- c(hap2.cis.concordance, hap2.trans.concordance)  
      #} else {
        mask.hap2 <- c(mask.hap2, i)
      }
    }

    ## Remove unreliable single cell haplotypes unless there is at least one left!!!
    if (length(mask.hap1) > 0 & length(mask.hap1) < nrow(hap1)) {
      hap1.files <- hap1.files[-mask.hap1]
      assembled.hap1 <- assembled.hap1[-mask.hap1]
      hap1 <- hap1[-mask.hap1, drop=FALSE,]
      hap1.quals <- hap1.quals[-mask.hap1, drop=FALSE,]
    }  
  
    if (length(mask.hap2) > 0 & length(mask.hap2) < nrow(hap2)) {
      hap2.files <- hap2.files[-mask.hap2]
      assembled.hap2 <- assembled.hap2[-mask.hap2]
      hap2 <- hap2[-mask.hap2, drop=FALSE,]
      hap2.quals <- hap2.quals[-mask.hap2, drop=FALSE,]
    }
  
    if (length(mask.hap1) == nrow(hap1) | length(mask.hap2) == nrow(hap2)) {
      warning("    Haplotypes did not pass concordance threshold: ", concordance, " !!!")
    }  
    
    ## Recalculate consensus haplotyes after filtering unreliable single cell haplotypes (using user defined settings)
    hap1.cons <- exportConsensus(data.bases = hap1, data.quals = hap1.quals, min.cov=1, translateBases=translateBases)
    hap2.cons <- exportConsensus(data.bases = hap2, data.quals = hap2.quals, min.cov=1, translateBases=translateBases)
  
    hap1.GenomicPos <- genomic.pos[hap1.cons$pos]
    hap2.GenomicPos <- genomic.pos[hap2.cons$pos]
  
    hap1.cons$pos <- hap1.GenomicPos
    hap2.cons$pos <- hap2.GenomicPos
  
  } else {
    
    assembled.hap1 <- list()
    assembled.hap2 <- list()
    assembled.haps <- list()
    for (i in 1:length(hap1.files)) {
      hap1.file <- hap1.files[i]
      hap2.file <- hap2.files[i]
      assembled.hap1[[hap1.file]] <- cbind(0, 0)
      assembled.hap2[[hap2.file]] <- cbind(0, 0)
    }  
    
    hap1.cons <- exportConsensus(hap1, hap1.quals, min.cov=1, translateBases=translateBases)
    hap2.cons <- exportConsensus(hap2, hap2.quals, min.cov=1, translateBases=translateBases)
    
    hap1.GenomicPos <- genomic.pos[hap1.cons$pos]
    hap2.GenomicPos <- genomic.pos[hap2.cons$pos]
    
    hap1.cons$pos <- hap1.GenomicPos
    hap2.cons$pos <- hap2.GenomicPos
    
    for (i in 1:length(hap1.files)) {
      ## export all phased regions (unfiltered)
      phased.info <- unlist(strsplit(hap1.files[i], '__'))
      ID <- phased.info[1]	
      hap1.phase <- phased.info[3]	
      phased.region <-  phased.info[2]
      phased.region <- unlist(strsplit(phased.region, ':|-'))		
      hap2.phase <- unlist(strsplit(hap2.files[i], '__'))[3]
      phase <- paste0(hap1.phase, hap2.phase)
      assembled.haps[[i]] <- paste(sep = "\t", c("NA", ID, phased.region, phase, 0, 0, 0, 0))
    }
    #assembled.haps[[i]] <- paste(sep = "\t", rep("NA", 8))		
  }  
  
  assembled.haps <- do.call(rbind, assembled.haps)	  

  ## Export data as a single data object (List)
  assem.haps <- list()
  assem.haps[['hap1.cons']] <- hap1.cons
  assem.haps[['hap2.cons']] <- hap2.cons
  assem.haps[['hap1.files']] <- assembled.hap1
  assem.haps[['hap2.files']] <- assembled.hap2
  assem.haps[['assem.haps']] <- assembled.haps
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  
  return(assem.haps)
}    
