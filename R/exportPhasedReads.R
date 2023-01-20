#' Export Strand-seq phased reads in FASTA or FASTQ format.
#'
#' This function will take phase information for each single cell and each haplotype informative region 
#' and will divide all mapped reads in BAM format into haplotype specific sets. In case of paired-end 
#' reads, read1 and read2 will be exported in separate files.
#' 
#' @param phased.tab A data table containing haplotype assignment for each single cell and haplotype informative region.
#' This table is standard output of \pkg{StrandPhaseR} and cen be found in subfolder 'Phased/phased_haps.txt'. 
#' @param by.matepair If set to \code{TRUE} (the default) paired-end reads will be exported per mate-pair (read1 and read2). 
#' @param format User defined output format, either "fasta" or "fastq" (the default). 
#' @inheritParams bamregion2GRanges
#' @inheritParams phaseChromosome
#' @importFrom Rsamtools indexBam ScanBamParam
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky
#' @export
#' 
exportPhasedReads <- function(phased.tab=NULL, inputfolder=NULL, outputfolder=NULL, by.matepair=TRUE, min.mapq=10, pairedEndReads=TRUE, format='fastq') {
  ## Set outputfolder
  if (!is.null(outputfolder)) {
    out.dir <- file.path(outputfolder, 'Phased_reads')
  } else {
    out.dir <- file.path(getwd(), 'Phased_reads')
  }
  ## Create output directory if not existing
  if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = TRUE)
  }
  
  ## Read in phasing table
  phased.df <- read.table(phased.tab, header = TRUE, stringsAsFactors = FALSE)
  phased.gr <- makeGRangesFromDataFrame(phased.df, keep.extra.columns = TRUE)
  
  ## Loop over all BAMs and export reads
  bamfiles <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
  for (i in seq_along(bamfiles)) {
    bamfile <- bamfiles[i]
    bam.name <- basename(bamfile)
    message("Processing bamfile: ", bam.name)
    ## Read in bam file ##
    ## Check if bamindex exists
    bamindex <- bamfile
    bamindex.raw <- sub('\\.bai$', '', bamindex)
    bamindex <- paste0(bamindex.raw,'.bai')
    if (!file.exists(bamindex)) {
      bamindex.own <- Rsamtools::indexBam(bamfile)
      warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
      bamindex <- bamindex.own
    }
    
    ## Read in raw reads
    bam.regions <- phased.gr[phased.gr$cell == bam.name]
    if (pairedEndReads) {
      suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(bam.regions), what=c('seq', 'qual','mapq','qname','flag'), flag=scanBamFlag(isProperPair=TRUE))) )	
    } else {
      suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(bam.regions), what=c('seq', 'qual','mapq','qname'), flag=scanBamFlag(isDuplicate=FALSE))) )
    } 

    ## Loop over chromosomal regions
    for (j in seq_along(bam.regions)) {
      bam.region <- bam.regions[j]
      data.raw.region <- subsetByOverlaps(data.raw, bam.region)
      chr.id <- as.character(seqnames(bam.region))
      ## Second mate of the pair will inherit directionality from the first mate of the pair
      if (pairedEndReads) {
        data.first <- as(GenomicAlignments::first(data.raw.region), 'GRanges')
        data.last <- as(GenomicAlignments::last(data.raw.region), 'GRanges')
        strand(data.last) <- strand(data.first)
        frags <- GenomicRanges::sort(c(data.first, data.last), ignore.strand=TRUE)
      } else {
        frags <- as(data.raw.region, 'GRanges')
      }
      ## Filter reads by mapping quality
      if (!is.null(min.mapq) & min.mapq > 0) {
        frags <- frags[mcols(frags)$mapq >= min.mapq]
        if (length(frags) == 0) {
          stop("No reads left to process after 'min.mapq' filtering, exiting!!!")
        }
      }
      
      ## Separate reads per haplotype
      if (bam.region$class == 'CW') {
        frags.h1 <- frags[strand(frags) == '+']
        frags.h2 <- frags[strand(frags) == '-']
      } else {
        frags.h1 <- frags[strand(frags) == '-']
        frags.h2 <- frags[strand(frags) == '+']
      }
      
      ## Export reads
      if (pairedEndReads) {
        ## Export pairedEndReads
        if (by.matepair) {
          ## Split reads by mate-pairs (read1 and read2)
          bit.flag <- bitwAnd(64, frags.h1$flag) ## Check first in pair
          r1.h1 <- frags.h1[bit.flag == 64]
          r2.h1 <- frags.h1[bit.flag == 0]
          bit.flag <- bitwAnd(64, frags.h2$flag)
          r1.h2 <- frags.h2[bit.flag == 64]
          r2.h2 <- frags.h2[bit.flag == 0]
          ## Read1 haplotype 1
          destination <- file.path(out.dir, paste0(chr.id, '_r1_h1.fastq.gz'))
          r1.h1.set <- r1.h1$seq
          names(r1.h1.set) <- r1.h1$qname
          Biostrings::writeXStringSet(x = r1.h1.set, append = TRUE, compress = TRUE, filepath = destination, qualities=r1.h1$qual, format = format)
          ## Read1 haplotype 2
          destination <- file.path(out.dir, paste0(chr.id, '_r1_h2.fastq.gz'))
          r1.h2.set <- r1.h2$seq
          names(r1.h2.set) <- r1.h2$qname
          Biostrings::writeXStringSet(x = r1.h2.set, append = TRUE, compress = TRUE, filepath = destination, qualities=r1.h2$qual, format = format)
          ## Read2 haplotype 1
          destination <- file.path(out.dir, paste0(chr.id, '_r2_h1.fastq.gz'))
          r2.h1.set <- r2.h1$seq
          names(r2.h1.set) <- r2.h1$qname
          Biostrings::writeXStringSet(x = r2.h1.set, append = TRUE, compress = TRUE, filepath = destination, qualities=r2.h1$qual, format = format)
          ## Read2 haplotype 2
          destination <- file.path(out.dir, paste0(chr.id, '_r2_h2.fastq.gz'))
          r2.h2.set <- r2.h2$seq
          names(r2.h2.set) <- r2.h2$qname
          Biostrings::writeXStringSet(x = r2.h2.set, append = TRUE, compress = TRUE, filepath = destination, qualities=r2.h2$qual, format = format)
        } else {
          ## Read1 haplotype 1
          destination <- file.path(out.dir, paste0(chr.id, '_r0_h1.fastq.gz'))
          frags.h1.set <- frags.h1$seq
          names(frags.h1.set) <- frags.h1$qname
          Biostrings::writeXStringSet(x = frags.h1.set, append = TRUE, compress = TRUE, filepath = destination, qualities=frags.h1$qual, format = format)
          ## Read1 haplotype 2
          destination <- file.path(out.dir, paste0(chr.id, '_r0_h2.fastq.gz'))
          frags.h2.set <- frags.h2$seq
          names(frags.h2.set) <- frags.h2$qname
          Biostrings::writeXStringSet(x = frags.h2.set, append = TRUE, compress = TRUE, filepath = destination, qualities=frags.h2$qual, format = format)
        }  
      } else {
        ## Export singleEndReads
        ## Read1 haplotype 1
        destination <- file.path(out.dir, paste0(chr.id, '_r0_h1.fastq.gz'))
        frags.h1.set <- frags.h1$seq
        names(frags.h1.set) <- frags.h1$qname
        Biostrings::writeXStringSet(x = frags.h1.set, append = TRUE, compress = TRUE, filepath = destination, qualities=frags.h1$qual, format = format)
        ## Read1 haplotype 2
        destination <- file.path(out.dir, paste0(chr.id, '_r0_h2.fastq.gz'))
        frags.h2.set <- frags.h2$seq
        names(frags.h2.set) <- frags.h2$qname
        Biostrings::writeXStringSet(x = frags.h2.set, append = TRUE, compress = TRUE, filepath = destination, qualities=frags.h2$qual, format = format)
      }
    }
  }
}
