#' This function calculates score of a matrix as a sum of partial scores of each columns
#' 
#' @param matrix object to calculate score for
#' @param covered.pos vector of snv positions to consider
#' 
#' @author David Porubsky
#' @export


calcMatrixScore <- function(matrix, covered.pos) {
  sub.m <- matrix[,covered.pos]
  sub.m <- sub.m[,colSums(sub.m) > 0]
  base.freq <- apply(sub.m, 2, function(x) table(x[x>0]))
  scores <- sapply(base.freq, function(x) sum(x)-max(x))
  return(sum(scores))
}



#' @author David Porubsky
#' @export

calcProb <- function(bases, quals) {
  probs <- vector()  
  for (i in 1:4) { #calculate phred score for each base (1:4)
    equal <- (1-10**(-quals[bases == i]/10)) # if base in bases if equals to tested base
    non.equal <- (1/3) * (10**(-quals[!bases == i]/10)) # if base in bases if not equals to tested base
    prob <- prod(equal, non.equal) # multiply calculated probabilities
    probs[i] <- prob    
  }
  probs <- probs/sum(probs) # normalize calculated probabilities by the sum of all probabilities
  cons <- which.max(probs) # get the base with max probability
  score <- probs[which.max(probs)]#  get the values of max probability
  return(list(cons,score))
}



#' @author David Porubsky
#' @export

calcEnt <- function(bases) {
  
  #log2 function (modif: return zero for zero values)
  getlog2 <- function(x) {
    if (x == 0 ) { 
      return(0)
    } else {  
      return(log2(x))
    }  
  }
  
  counts <- c(0,0,0,0)
  
  # count frequency of each base in a column 
  for (i in 1:4) {
    counts[i] <- length(bases[bases == i])
  }
  probs <- counts/sum(counts) # get probabilities
  
  ent <- -( probs[1]*getlog2(probs[1]) + probs[2]*getlog2(probs[2]) + probs[3]*getlog2(probs[3]) + probs[4]*getlog2(probs[4]) )
  return(ent)    
}