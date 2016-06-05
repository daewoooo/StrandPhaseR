#' This will take phased info for each haplotype and will split directional reads of each single cell into haplotype
#' specific reads
#' 
#' @param data.object containing sorted and filtered watson and crick haplotypes of each single cell along with phase info
#' @param gr.list
#' 
#' @author David Porubsky
#' @export

splitReads <- function(data.object=NULL, bamfilespath=bamfilespath, pairedEndReads=FALSE, min.mapq=10) {
  message("Printing haplotypes")
  ptm <- proc.time()
  
  hap1.files <- data.object[['hap1.files']]
  hap2.files <- data.object[['hap2.files']]
  
  hap1 <- names(hap1.files)
  hap2 <- names(hap2.files)
  
  ## reformat phased file lists for both haplotypes
  hap1.parts <- unlist(strsplit(hap1, "__"))
  filenames.hap1 <- hap1.parts[seq(1, length(hap1.parts), 3)]
  regions <- hap1.parts[seq(2, length(hap1.parts), 3)]
  phase <- hap1.parts[seq(3, length(hap1.parts), 3)]
  hap1.list <- setNames(as.list(phase), paste(filenames.hap1, regions, sep="__"))
  
  #hap1.list <- mapply(c, regions, phase, SIMPLIFY = FALSE)
  #names(hap1.list) <- filenames.hap1
  
  hap2.parts <- unlist(strsplit(hap2, "__"))
  filenames.hap2 <- hap2.parts[seq(1, length(hap2.parts), 3)]
  regions <- hap2.parts[seq(2, length(hap2.parts), 3)]
  phase <- hap2.parts[seq(3, length(hap2.parts), 3)]
  hap2.list <- setNames(as.list(phase), paste(filenames.hap2, regions, sep="__"))
  
  #hap2.list <- mapply(c, regions, phase, SIMPLIFY = FALSE)
  #names(hap2.list) <- filenames.hap2
  
  region.IDs <- union(names(hap1.list), names(hap2.list))
  
  hap1 <- GRangesList()
  hap2 <- GRangesList()
  for (region in region.IDs) {
   
    ##split region identifier
    parts <- unlist(strsplit(region, "__"))
    filename <- parts[1]
    parts2 <- unlist(strsplit(parts[2], ":"))
    chr <- parts2[1]
    parts3 <- unlist(strsplit(parts2[2], "-"))
    start <- parts3[1]
    end <- parts3[2]
    
    region.gr <- GRanges(chr, IRanges(start=as.numeric(start), end=as.numeric(end)))
    
    #get data from gr.list
    #data <- gr.list[[filename]]
    #data <- data[,3] #take along mapq metadata column from mcols 1,2,3,4
    bamfile <- file.path(bamfilespath, filename)
    data <- bamregion2GRanges(bamfile, region=region.gr, pairedEndReads=pairedEndReads, min.mapq=min.mapq)
    data <- data[,3] #take along mapq metadata column from mcols 1,2,3,4
    
    ## Split reads in WC regions belonging to the Hap1
    if (exists(region, where=hap1.list)) {
      phase <- hap1.list[[region]]
      if (phase == 'W') {
        minus.reads <- data[strand(data) == "-"]
        hap1[[region]] <- minus.reads
      } else {
        plus.reads <- data[strand(data) == "+"]
        hap1[[region]] <- plus.reads
      }
    }
    
    ## Split reads in WC regions belonging to the Hap2
    if (exists(region, where=hap2.list)) {
      phase <- hap2.list[[region]]
      if (phase == 'W') {
        minus.reads <- data[strand(data) == "-"]
        hap2[[region]] <- minus.reads
      } else {
        plus.reads <- data[strand(data) == "+"]
        hap2[[region]] <- plus.reads
      }
    }
  }
  
  ## Store data in GRanges object
  haps <- GenomicRanges::GRangesList()
  hap1 <- GenomicRanges::sort(unname(GenomicRanges::unlist(hap1)))
  hap2 <- GenomicRanges::sort(unname(GenomicRanges::unlist(hap2)))
  haps[['hap1']] <- hap1
  haps[['hap2']] <- hap2
  
  time <- proc.time() - ptm
  message("Time spent: ",round(time[3],2),"s")
  
  return(haps)
}  
