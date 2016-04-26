#' This funcion will load initial watson and crick matrices and will sort them according to the phase information
#' of each single cell
#' 
#' @param data.object containing watson and crick haplotypes of each single cell
#' @param num.iterations number of iteration to sort watson and crick matrices
#' @importFrom Biostrings alphabetFrequency
#' 
#' @author David Porubsky
#' @export


sortMatrices <- function(data.object, num.iterations=2) {
  crick.m <- data.object[['crick.bases']]
  watson.m <- data.object[['watson.bases']]
  crickQuals.m <- data.object[['crick.quals']]
  watsonQuals.m <- data.object[['watson.quals']]
  filename.IDs <- data.object[['row.IDs']]
  genomic.pos <- data.object[['genomic.pos']]
  
  files.hap1 <- paste(filename.IDs, "C", sep="__")
  files.hap2 <- paste(filename.IDs, "W", sep="__")
  
  ## store phased files after each iteration
  hap1 <- as.list(rep("", length(filename.IDs)))
  hap2 <- as.list(rep("", length(filename.IDs)))
  
  ptm <- proc.time()
  for (iter in 1:num.iterations) {
    message("Sorting matrices: iteration ", iter)
    for ( i in 1:length(filename.IDs)) {
    
      crick.pos <- which(crick.m[i,] > 0)
      watson.pos <- which(watson.m[i,] > 0)
      cov.pos <- union(crick.pos, watson.pos)
    
      ## initialize score
      crick.m.score <- calcMatrixScore(crick.m, cov.pos)
      watson.m.score <- calcMatrixScore(watson.m, cov.pos)
    
      ## swap rows in matrices
      hap1.filename <- files.hap1[i]
      hap2.filename <- files.hap2[i]
      files.hap1[i] <- hap2.filename
      files.hap2[i] <- hap1.filename
      
      hap1[[i]][iter] <- hap2.filename
      hap2[[i]][iter] <- hap1.filename
    
      crick.row <- crick.m[i,]
      watson.row <- watson.m[i,]
      crick.m[i,] <- watson.row
      watson.m[i,] <- crick.row
      
      crickQuals.row <- crickQuals.m[i,]
      watsonQuals.row <- watsonQuals.m[i,]
      crickQuals.m[i,] <- watsonQuals.row
      watsonQuals.m[i,] <- crickQuals.row
    
      ## calculate new score
      curr.crick.m.score <- calcMatrixScore(crick.m, cov.pos)
      curr.watson.m.score <- calcMatrixScore(watson.m, cov.pos)
    
      ## compare previous matrix score with score after swapping rows
      if ( (crick.m.score + watson.m.score) < (curr.crick.m.score + curr.watson.m.score) ) {
        files.hap1[i] <- hap1.filename
        files.hap2[i] <- hap2.filename
        
        crick.m[i,] <- crick.row
        watson.m[i,] <- watson.row
        crickQuals.m[i,] <- crickQuals.row
        watsonQuals.m[i,] <- watsonQuals.row
        
        hap1[[i]][iter] <- hap1.filename
        hap2[[i]][iter] <- hap2.filename
      }
    }  
  }
  sorted.matrices <- list()
  sorted.matrices[['hap1.bases']] <- crick.m
  sorted.matrices[['hap2.bases']] <- watson.m
  sorted.matrices[['hap1.quals']] <- crickQuals.m
  sorted.matrices[['hap2.quals']] <- watsonQuals.m
  sorted.matrices[['hap1.files']] <- files.hap1
  sorted.matrices[['hap2.files']] <- files.hap2
  sorted.matrices[['genomic.pos']] <- genomic.pos
  
  time <- proc.time() - ptm
  message("Time spent: ",round(time[3],2),"s")
  
  return(sorted.matrices)
}  
  
  
 