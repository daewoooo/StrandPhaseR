#' This function will take both assembled haplotypes and will try to fill gaps at positions where only one allele is phased
#' Such position have to be heterozygous so the alternative allele at this position can be reliably distinguished
#' 
#' @param data.object ...
#' @param merged.bam ...
#' @param chunkSize ...
#' @inheritParams bamregion2GRanges
#' @inheritParams phaseChromosome
#' @import GenomicRanges
#' @importFrom GenomicAlignments pileLettersAt
#' @importFrom Biostrings DNAStringSet
#' 
#' @author David Porubsky
#' @export

fillGapsWithBam <- function(data.object, merged.bam, min.mapq=10, min.baseq=30, translateBases=FALSE, chromosome=NULL, chunkSize=10000000, filterAltAlign=TRUE) {
 
  message(" Filling gaps in haplotypes", appendLF=FALSE); ptm <- proc.time()
  
  hap1.cons <- data.object[['hap1.cons']]
  hap2.cons <- data.object[['hap2.cons']]
  assembled.hap1 <- data.object[['hap1.files']]
  assembled.hap2 <- data.object[['hap2.files']]
  haps <- data.object[['assem.haps']]
  
  cov.both <- intersect(hap1.cons$pos,hap2.cons$pos) #get positions covered on both haplotypes
  
  hap1.gaps <- hap2.cons$pos[!hap2.cons$pos %in% cov.both] #postions covered only on hap2 = gaps in hap1
  hap2.gaps <- hap1.cons$pos[!hap1.cons$pos %in% cov.both] #postions covered only on hap1 = gaps in hap2

  hap1Gaps.gr <- GRanges(seqnames=as.character(chromosome), IRanges(start=hap1.gaps, end=hap1.gaps))
  hap2Gaps.gr <- GRanges(seqnames=as.character(chromosome), IRanges(start=hap2.gaps, end=hap2.gaps))  
  
  #initialize list to store HAP1 data
  hap1fillBases <- as.list(rep("", length(hap1.gaps)))
  hap1fillCov <- as.list(rep(0, length(hap1.gaps)))
  hap1fillQuals <- as.list(rep(0, length(hap1.gaps)))
  hap1fillscore <- as.list(rep(0, length(hap1.gaps)))
  hap1fillentropy <- as.list(rep(0, length(hap1.gaps)))
  
  #initialize list to store HAP1 data
  hap2fillBases <- as.list(rep("", length(hap2.gaps)))
  hap2fillCov <- as.list(rep(0, length(hap2.gaps)))
  hap2fillQuals <- as.list(rep(0, length(hap2.gaps)))
  hap2fillscore <- as.list(rep(0, length(hap2.gaps)))
  hap2fillentropy <- as.list(rep(0, length(hap2.gaps)))
  
  ## Get chromosome length
  file.header <- Rsamtools::scanBamHeader(merged.bam)[[1]]
  chrom.lengths <- file.header$targets
  chrom.len <- chrom.lengths[chromosome]
  #chrom.range <- GenomicRanges::GRanges(seqnames=chromosome, ranges=IRanges(start=1, end=chrom.len))
  
  #process data in chunks to save memory
  chunks <- unlist(tileGenome(chrom.len, tilewidth = chunkSize))

  for (i in 1:length(chunks)) {
    chunk <- chunks[i]

    data <- bamregion2GRanges(bamfile=merged.bam, bamindex=merged.bam, region=chunk, pairedEndReads=FALSE, min.mapq=min.mapq, filterAltAlign=filterAltAlign)
  
    ## extract read sequences
    pile.seq <- mcols(data)$seq

    ## extract base qualities
    pile.qual <- mcols(data)$qual

    hap1Gaps.piles <- GenomicAlignments::pileLettersAt(pile.seq, seqnames(data), start(data), mcols(data)$cigar, hap1Gaps.gr)
    hap2Gaps.piles <- GenomicAlignments::pileLettersAt(pile.seq, seqnames(data), start(data), mcols(data)$cigar, hap2Gaps.gr)

    hap1Gaps.quals <- GenomicAlignments::pileLettersAt(pile.qual, seqnames(data), start(data), mcols(data)$cigar, hap1Gaps.gr)
    hap2Gaps.quals <- GenomicAlignments::pileLettersAt(pile.qual, seqnames(data), start(data), mcols(data)$cigar, hap2Gaps.gr)

    ## Filter bases based on base quality
    df.hap1GapsPiles <- as(hap1Gaps.piles, "data.frame")
    df.hap2GapsPiles <- as(hap2Gaps.piles, "data.frame")
    df.hap1GapsQuals <- as(hap1Gaps.quals, "data.frame")
    df.hap2GapsQuals <- as(hap2Gaps.quals, "data.frame")

    bases.hap1Gaps <- strsplit(df.hap1GapsPiles$x, "")
    bases.hap2Gaps <- strsplit(df.hap2GapsPiles$x, "")
  
    ## translate raw base qualities into a number (for Sanger 33)
    quals.hap1Gaps <- sapply(df.hap1GapsQuals$x, function(x) as.numeric(charToRaw(x))-33)
    quals.hap2Gaps <- sapply(df.hap2GapsQuals$x, function(x) as.numeric(charToRaw(x))-33)
  
    filtbases.hap1Gaps <- mapply(function(X,Y) { X[Y >= min.baseq] }, X=bases.hap1Gaps, Y=quals.hap1Gaps)
    filtquals.hap1Gaps <- sapply(quals.hap1Gaps, function(x) x[x >= min.baseq])
    filtbases.hap2Gaps <- mapply(function(X,Y) { X[Y >= min.baseq] }, X=bases.hap2Gaps, Y=quals.hap2Gaps)
    filtquals.hap2Gaps <- sapply(quals.hap2Gaps, function(x) x[x >= min.baseq])
  
    #get bases phased on opposing haplotype
    knownBases.hap1Gaps <- hap2.cons[hap2.cons$pos %in% hap1.gaps,]$bases
    knownBases.hap2Gaps <- hap1.cons[hap1.cons$pos %in% hap2.gaps,]$bases
  
    ## Fill gaps for haplotype 1
    for (i in 1:length(filtbases.hap1Gaps)) { #loop over bases prefiltered based on base quality
      base2fill <- filtbases.hap1Gaps[[i]] #get bases
      qual2fill <- quals.hap1Gaps[[i]] #get base qualities
      phasedBase <- knownBases.hap1Gaps[i] #get phased base for this position
      bases <- base2fill[base2fill != phasedBase] #filter out bases coresponding to phased bases
      quals <- qual2fill[base2fill != phasedBase] #filter out base qualities coresponding to phased bases
      if (length(bases)) { #if position heterozygous => alternative base than phased base is present
        base.freq <- table(bases) #calc base frequency
      
        if (length(base.freq) < 2) { #only one alternative base present
          max.base <- names(which.max(base.freq)) #get max abundant base
          quals <- quals[bases == max.base] #get base qualities for max.base
          bases <- bases[bases == max.base] #get bases equal to max.base
          bases <- chartr("ACGT", "1234", bases) #transfer bases to numbers
          score <- calcProb(bases, quals) #calculate probability score
          hap1fillBases[[i]] <- max.base
          hap1fillCov[[i]] <- length(bases)  
          hap1fillQuals[[i]] <- round(mean(quals))
          hap1fillscore[[i]] <- score[[2]]
        } else { #more than one alternative base prensent
          if (max(base.freq) > min(base.freq)) { #one alternative base has to be more abundant
            max.base <- names(which.max(base.freq))
            quals <- quals[bases == max.base]
            bases <- bases[bases == max.base]
            bases <- chartr("ACGT", "1234", bases)
            score <- calcProb(bases, quals)
            hap1fillBases[[i]] <- max.base
            hap1fillCov[[i]] <- length(bases)
            hap1fillQuals[[i]] <- round(mean(quals))
            hap1fillscore[[i]] <- score[[2]]
          }  
        }
      }
    }
  
    ## Fill gaps for haplotype 2
    for (i in 1:length(filtbases.hap2Gaps)) {
      base2fill <- filtbases.hap2Gaps[[i]]
      qual2fill <- quals.hap2Gaps[[i]]
      phasedBase <- knownBases.hap2Gaps[i]
      bases <- base2fill[base2fill != phasedBase]
      quals <- qual2fill[base2fill != phasedBase]
      if (length(bases)) {
        base.freq <- table(bases)
      
        if (length(base.freq) < 2) {
          max.base <- names(which.max(base.freq))
          quals <- quals[bases == max.base]
          bases <- bases[bases == max.base]
          bases <- chartr("ACGT", "1234", bases)
          score <- calcProb(bases, quals)
          hap2fillBases[[i]] <- max.base
          hap2fillCov[[i]] <- length(bases)  
          hap2fillQuals[[i]] <- round(mean(quals))
          hap2fillscore[[i]] <- score[[2]]
        } else {
          if (max(base.freq) > min(base.freq)) {
            max.base <- names(which.max(base.freq))
            quals <- quals[bases == max.base]
            bases <- bases[bases == max.base]
            bases <- chartr("ACGT", "1234", bases)
            score <- calcProb(bases, quals)
            hap2fillBases[[i]] <- max.base
            hap2fillCov[[i]] <- length(bases)
            hap2fillQuals[[i]] <- round(mean(quals))
            hap2fillscore[[i]] <- score[[2]]
          }  
        }
      }
    }
  } 

  #Summarize filled snvs for Hap1    
  mask <- hap1fillBases != ''
  filledHap1.bases <- unlist( hap1fillBases[mask] )
  filledHap1.cov <- unlist( hap1fillCov[mask] )
  filledHap1.quals <- unlist( hap1fillQuals[mask] )
  filledHap1.score <- unlist( hap1fillscore[mask] )
  filledHap1.entropy <- unlist( hap1fillentropy[mask] )
  filledHap1.pos <- hap1.gaps[mask]
  
## translate calcualted probabilities of correcty base call to base qualities
#if (score2qual) {
#  prob.err <- filledHap1.score - 1 #get probability of error base call (score=probability of correct base call)
#  prob.err[prob.err < 0.00006] <- 0.00006 #do this if the error probability is lower than min possible
#  filledHap1.score <- round( log10(prob.err)*-10 )
#}
  
  if (translateBases) {
    filledHap1.bases <- chartr("1234", "ACGT", filledHap1.bases)
  }
  
  filled.hap1 <- data.frame(pos=filledHap1.pos, bases=filledHap1.bases, cov=filledHap1.cov, score=filledHap1.score, ent=filledHap1.entropy)
  filled.hap1 <- rbind(hap1.cons, filled.hap1)
  filled.hap1 <- filled.hap1[order(filled.hap1$pos),]
  rownames(filled.hap1) <- NULL
  
  #Summarize filled snvs for Hap2    
  mask <- hap2fillBases != ''
  filledHap2.bases <- unlist( hap2fillBases[mask] )
  filledHap2.cov <- unlist( hap2fillCov[mask] )
  filledHap2.quals <- unlist( hap2fillQuals[mask] )
  filledHap2.score <- unlist( hap2fillscore[mask] )
  filledHap2.entropy <- unlist( hap2fillentropy[mask] )
  filledHap2.pos <- hap2.gaps[mask]
  
## translate calcualted probabilities of correcty base call to base qualities
#if (score2qual) {
#  prob.err <- filledHap2.score - 1 #get probability of error base call (score=probability of correct base call)
#  prob.err[prob.err < 0.00006] <- 0.00006 #do this if the error probability is lower than min possible
#  filledHap2.score <- round( log10(prob.err)*-10 )
#}
  
  if (translateBases) {
    filledHap2.bases <- chartr("1234", "ACGT", filledHap2.bases)
  }  
  
  filled.hap2 <- data.frame(pos=filledHap2.pos, bases=filledHap2.bases, cov=filledHap2.cov, score=filledHap2.score, ent=filledHap2.entropy)
  filled.hap2 <- rbind(hap2.cons, filled.hap2)
  filled.hap2 <- filled.hap2[order(filled.hap2$pos),]
  rownames(filled.hap2) <- NULL
  
  ## Export data as a single data object (List)
  filled.haps <- list()
  filled.haps[['hap1.cons']] <- filled.hap1
  filled.haps[['hap2.cons']] <- filled.hap2
  filled.haps[['hap1.files']] <- assembled.hap1
  filled.haps[['hap2.files']] <- assembled.hap2
  filled.haps[['assem.haps']] <- haps
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  return(filled.haps)
}


#' This function will take both assembled haplotypes and will try to fill gaps at positions where only one allele is phased
#' Such position have to be heterozygous so the alternative allele at this position can be reliably distinguished
#' 
#' @param data.object ...
#' @param ref.vcf ...
#' @inheritParams phaseChromosome
#' 
#' @author David Porubsky
#' @export

fillGapsWithVCF <- function(data.object, ref.vcf, chromosome=NULL) {
  
  message(" Filling gaps in haplotypes", appendLF=FALSE); ptm <- proc.time()
  
  hap1.cons <- data.object[['hap1.cons']]
  hap2.cons <- data.object[['hap2.cons']]
  assembled.hap1 <- data.object[['hap1.files']]
  assembled.hap2 <- data.object[['hap2.files']]
  haps <- data.object[['assem.haps']]
  
  cov.both <- intersect(hap1.cons$pos,hap2.cons$pos) #get positions covered on both haplotypes
  
  hap1.gaps <- hap2.cons$pos[!hap2.cons$pos %in% cov.both] #postions covered only on hap2 = gaps in hap1
  hap2.gaps <- hap1.cons$pos[!hap1.cons$pos %in% cov.both] #postions covered only on hap1 = gaps in hap2
  hap1Gaps.gr <- GRanges(seqnames=as.character(chromosome), IRanges(start=hap1.gaps, end=hap1.gaps))
  hap2Gaps.gr <- GRanges(seqnames=as.character(chromosome), IRanges(start=hap2.gaps, end=hap2.gaps))  
  
  #read in referene vcf file
  suppressMessages( snvs <- vcf2ranges(vcfFile=ref.vcf, genotypeField=1, chromosome=chromosome) )
  
  ## Remove duplicated SNV positions if exists
  msg <- NULL
  dup.snvs <- any(duplicated(snvs))
  if (dup.snvs) {
    mask <- !duplicated(snvs)
    removed.snvs <- table(mask)['FALSE']
    snvs <- snvs[mask]
    msg <- paste0("    Removed ", removed.snvs, " duplicated SNV positions!!!")
  }
  
  # Fill gaps in haplotype 1
  hits.hap1 <- IRanges::findOverlaps(snvs, hap1Gaps.gr)
  missed.hap1 <- snvs[S4Vectors::queryHits(hits.hap1)]
  known.hap2 <- as.character(hap2.cons$bases[match(start(missed.hap1), hap2.cons$pos)])
  missed.hap1$fill <- rep("", length(missed.hap1))
  missed.hap1$fill[known.hap2 == missed.hap1$ref] <- missed.hap1$alt[known.hap2 == missed.hap1$ref]
  missed.hap1$fill[known.hap2 == missed.hap1$alt] <- missed.hap1$ref[known.hap2 == missed.hap1$alt]
  
  filled.hap1 <- data.frame(pos=as.numeric(start(missed.hap1)), bases=missed.hap1$fill, cov=0, score=0, ent=0)
  filled.hap1 <- filled.hap1[filled.hap1$bases != '',] #filter unfilled positions
  filled.hap1 <- rbind(hap1.cons, filled.hap1)
  filled.hap1 <- filled.hap1[order(filled.hap1$pos),]
  rownames(filled.hap1) <- NULL
  
  # Fill gaps in haplotype 2
  hits.hap2 <- findOverlaps(snvs, hap2Gaps.gr)
  missed.hap2 <- snvs[S4Vectors::queryHits(hits.hap2)]
  known.hap1 <- as.character(hap1.cons$bases[match(start(missed.hap2), hap1.cons$pos)])
  missed.hap2$fill <- rep("", length(missed.hap2))
  missed.hap2$fill[known.hap1 == missed.hap2$ref] <- missed.hap2$alt[known.hap1 == missed.hap2$ref]
  missed.hap2$fill[known.hap1 == missed.hap2$alt] <- missed.hap2$ref[known.hap1 == missed.hap2$alt]
  
  filled.hap2 <- data.frame(pos=as.numeric(start(missed.hap2)), bases=missed.hap2$fill, cov=0, score=0, ent=0)
  filled.hap2 <- filled.hap2[filled.hap2$bases != '',] #filter unfilled positions
  filled.hap2 <- rbind(hap2.cons, filled.hap2)
  filled.hap2 <- filled.hap2[order(filled.hap2$pos),]
  rownames(filled.hap2) <- NULL
  
  ## Export data as a single data object (List)
  filled.haps <- list()
  filled.haps[['hap1.cons']] <- filled.hap1
  filled.haps[['hap2.cons']] <- filled.hap2
  filled.haps[['hap1.files']] <- assembled.hap1
  filled.haps[['hap2.files']] <- assembled.hap2
  filled.haps[['assem.haps']] <- haps
  
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  if (is.character(msg)) {
    message(msg)
  }
  return(filled.haps)
}
