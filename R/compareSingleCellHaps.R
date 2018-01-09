#' This will take sorted matrices and will calculate concordance of each single cell to the consensus haplotypes in order to assemble highly 
#' accurate haplotypes
#' 
#' @param consensusHaps A \code{list} returned by function \code{\link{assembleHaps}}
#' @param sortedHaps  A \code{list} returned by function \code{\link{sortMatrices}}
#' @param bin.size A size of the window measured in number of SNV to scan single cell haplotypes.
#' @param step A step (number of SNVs) to move the window.
#' @importFrom zoo rollapply
#' 
#' @author David Porubsky
#' @export

compareSingleCellHaps <- function(consensusHaps=NULL, sortedHaps=NULL, bin.size=10, step=1) {
  message(" Scanning single-cell haplotypes", appendLF=F); ptm <- proc.time()

  #helper function
  getSingleCellSimil <- function(df) {
    cell.hap <- df[,1]
    hap1.cons <- df[,2]
    hap2.cons <- df[,3]
    #cell.hap.string <- paste(chartr("1234", "ACGT", cell.hap), collapse = "") [SLOW]
    #hap1.cons.string <- paste(chartr("1234", "ACGT", hap1.cons), collapse = "") [SLOW]
    #hap2.cons.string <- paste(chartr("1234", "ACGT", hap2.cons), collapse = "") [SLOW]
    agree.hap1 <- length(cell.hap[cell.hap == hap1.cons])
    agree.hap1.perc <- as.numeric(agree.hap1/length(cell.hap))
    agree.hap2 <- length(cell.hap[cell.hap == hap2.cons])
    agree.hap2.perc <- as.numeric(agree.hap2/length(cell.hap))
    
    #df.new <- data.frame(start=min(as.numeric(rownames(df))), end=max(as.numeric(rownames(df))), cell.hap=cell.hap.string, hap1.cons=hap1.cons.string, hap2.cons=hap2.cons.string, cons1.simil=agree.hap1.perc, cons2.simil=agree.hap2.perc)
    df.new <- data.frame(start=min(as.numeric(rownames(df))), end=max(as.numeric(rownames(df))), cons1.simil=agree.hap1.perc, cons2.simil=agree.hap2.perc)
    return(df.new)
  }
    
  #get the required data
  hap1.cons <- consensusHaps$hap1.cons
  hap2.cons <- consensusHaps$hap2.cons
  hap1.cells.srt <- sortedHaps$hap1.bases 
  hap2.cells.srt <- sortedHaps$hap2.bases
  hap1.cells.names <- sortedHaps$hap1.files 
  hap2.cells.names <- sortedHaps$hap2.files
  gen.positions <- sortedHaps$genomic.pos
  
  #convert consensus bases to numbers
  hap1.cons.code <- as.numeric(chartr("ACGT", "1234", hap1.cons$bases))
  hap2.cons.code <- as.numeric(chartr("ACGT", "1234", hap2.cons$bases))
  names(hap1.cons.code) <- hap1.cons$pos
  names(hap2.cons.code) <- hap2.cons$pos
  
  #check if there is enough SNV variants to scan single cell haplotypes (at least 10 times more)
  if (nrow(hap1.cons)>10*bin.size & nrow(hap1.cons)>10*bin.size) {
  
    cell.comparisons <- list()
    cell.haps.allSNVs <- list()
    for (i in 1:nrow(hap1.cells.srt)) {
      cell.hap1 <- hap1.cells.srt[i,]
      cell.hap2 <- hap2.cells.srt[i,]
      names(cell.hap1) <- gen.positions
      names(cell.hap2) <- gen.positions
      cell.ID <- strsplit(hap1.cells.names[i], split = "\\.")[[1]][1]
    
      #filter only covered SNV in a given cell and given haplotype
      cell.hap1 <- cell.hap1[cell.hap1 != 0]
      cell.hap2 <- cell.hap2[cell.hap2 != 0]
      #compara haps per SNV position
      #cell.hap1.comp <- rep(0, length(cell.hap1))
      #cell.hap1.comp[ cell.hap1 == hap1.cons.code[names(cell.hap1)] ] <- 1
      #cell.hap1.comp[ cell.hap1 == hap2.cons.code[names(cell.hap1)] ] <- 2
      #cell.hap2.comp <- rep(0, length(cell.hap2))
      #cell.hap2.comp[ cell.hap2 == hap1.cons.code[names(cell.hap2)] ] <- 1
      #cell.hap2.comp[ cell.hap2 == hap2.cons.code[names(cell.hap2)] ] <- 2
    
      #prepara data frame structure to compare single cell haplotypes
      hap1.comp <- data.frame(cell.hap1=cell.hap1, cons.hap1=hap1.cons.code[names(cell.hap1)], cons.hap2=hap2.cons.code[names(cell.hap1)])
      hap1.comp <- hap1.comp[complete.cases(hap1.comp),]
      hap2.comp <- data.frame(cell.hap2=cell.hap2, cons.hap1=hap1.cons.code[names(cell.hap2)], cons.hap2=hap2.cons.code[names(cell.hap2)])
      hap2.comp <- hap2.comp[complete.cases(hap2.comp),]
    
      #per position haplotype comparison
      #hap1.perPos.comp <- re
      #hap1.comp$cell.hap1 == hap1.comp$cons.hap1
    
      #binned haplotype comparison
      if (nrow(hap1.comp) >= bin.size) {
        cons1.simil <- zoo::rollapply(hap1.comp, width = bin.size, by=step, by.column = F, FUN=getSingleCellSimil)
      } else {
        cons1.simil <- getSingleCellSimil(hap1.comp)
      }  
    
      if (nrow(hap2.comp) >= bin.size) {
        cons2.simil <- zoo::rollapply(hap2.comp, width = bin.size, by=step, by.column = F, FUN=getSingleCellSimil) 
      } else {
        cons2.simil <- getSingleCellSimil(hap2.comp)
      }  
    
      cons1.simil <- as.data.frame(cons1.simil) 
      cons2.simil <- as.data.frame(cons2.simil) 
      cons1.simil$hap <- 'H1'
      cons2.simil$hap <- 'H2'
      cell.comp <- rbind(cons1.simil, cons2.simil)
      cell.comp$CellID <- paste(cell.ID, cell.comp$hap, sep = "_")
      cell.comparisons[[i]] <- cell.comp
    }  
    
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
    return(cell.comparisons)
  
  } else {
    message("\n Insufficient SNV density to scan single-cell haplotypes, skipping ...")
    return(NULL)
  }    
  #cell.comparisons.df <- do.call(rbind, cell.comparisons)
}