#' This funcion will load sorted watson and crick matrices and will collapsed them in order to obtain consensus string of each matrix
#' 
#' @param data.object containing sorted watson and crick haplotypes of each single cell
#' @param min.cov dtermines minimal number of cells covering each position to consider in consensus
#' @param translateBases translates integer coded bases (1,2,3,4) into letters (A,C,G,T) 
#' 
#' @author David Porubsky
#' @export

exportConsensus <- function(data.bases, data.quals, filt.ambig=TRUE, min.cov=2, translateBases=FALSE, score2qual=FALSE) {
  #base.freq <- apply(data.bases, 2, function(x) table(x[x>0])) ## TODO avoid warnings here
  indices <- which(data.bases!=0,arr.ind = T) # get indices of non-zero values
  values <- data.bases[indices] # get non-zero values based on indices
  col.vals <- split(values, (indices[,2])) # split values based on column indices
  
  indices <- which(data.quals!=0,arr.ind = T) # get indices of non-zero values
  values <- data.quals[indices] # get non-zero values based on indices
  col.quals <- split(values, (indices[,2])) # split values based on column indices 
  
  base.freq <- lapply(col.vals, table) # get base frequencies for each column
  
  positions <- as.numeric(names(base.freq)) # get position of each snv/column
  
  ## Filter SNV positions (columns) which do not pass set criteria (coverge, ambiguity)
  if (filt.ambig) {
    max.cov <- sapply(base.freq, function(x) max(x)) # get coverage of highest covered base
    score <- sapply(base.freq, function(x) sum(x)-max(x)) 
    mask <- which(max.cov > score & max.cov >= min.cov) # consider only bases passing filtering criteria
    base.freq <- base.freq[names(mask)]
    col.quals <- col.quals[names(mask)]
    col.vals <- col.vals[names(mask)]
    positions <- as.numeric(names(mask))
  }
  
  ## Calculate entropy values for each column
  entropy <- c(rep(0, length(positions)))
  entropy <- sapply(col.vals, calcEnt)
  
  ## Calculate Phred score values for each column
  scores <- c(rep(0, length(positions)))

  phredScore <- mapply(calcProb, bases=col.vals, quals=col.quals)
  bases <- unlist(phredScore[1,])
  scores <- unlist(phredScore[2,])
  
  ## translate calcualted probabilities of correcty base call to base qualities
  if (score2qual) {
    prob.err <- scores - 1 #get probability of error base call (score=probability of correct base call)
    prob.err[prob.err < 0.00006] <- 0.00006 #do this if the error probability is lower than min possible
    scores <- round( log10(prob.err)*-10 )
  }
  
  if (translateBases) {
    bases <- chartr("1234", "ACGT", bases)
  }  
  
  cov = sapply(col.vals, length)
  
  assem.haps <- data.frame(pos=positions, bases=bases, cov=cov, score=scores, ent=entropy)
  rownames(assem.haps) <- NULL
  return(assem.haps)
}