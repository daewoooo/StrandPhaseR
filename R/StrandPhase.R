#' Wrapper function 
#'
#' This function will move through .bam files in a folder and perform several steps (see Details).
#'
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @inheritParams phaseChromosome
#' @import foreach
#' @import doParallel 

#' @author David Porubsky
#' @export

strandPhaseR <- function(inputfolder, outputfolder='./StrandPhaseR_analysis', configfile=NULL, numCPU=1, positions=NULL, WCregions=NULL, chromosomes=NULL, pairedEndReads=TRUE, min.mapq=10, min.baseq=30, num.iterations=2, translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, callBreaks=FALSE, exportVCF=NULL, bsGenome=NULL) {
  
  #=======================
  ### Helper functions ###
  #=======================
  as.object <- function(x) {
    return(eval(parse(text=x)))
  }
  
  #========================
  ### General variables ###
  #========================
  conf <- NULL
  if (is.character(configfile)) {
    ## Read config file ##
    errstring <- tryCatch({
      conf <- readConfig(configfile)
      errstring <- ''
    }, error = function(err) {
      errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
      stop(errstring)
    }
  }
  
  ## Convert BSgenome to string if necessary
  if (class(bsGenome)=='BSgenome') {
    bsGenome <- attributes(bsGenome)$pkgname
  }
  
  ## Put options into list and merge with conf
  params <- list(numCPU=numCPU, positions=positions, WCregions=WCregions, chromosomes=chromosomes, pairedEndReads=pairedEndReads, min.mapq=min.mapq, min.baseq=min.baseq, num.iterations=num.iterations, translateBases=translateBases, fillMissAllele=fillMissAllele, splitPhasedReads=splitPhasedReads, compareSingleCells=compareSingleCells, callBreaks=callBreaks, exportVCF=exportVCF, bsGenome=bsGenome)
  conf <- c(conf, params[setdiff(names(params),names(conf))])
  
  #===================
  ### Input checks ###
  #===================
  ## Check user input
  if (!is.null(conf[['exportVCF']]) & is.null(conf[['bsGenome']])) {
    warning("VCf file cannot be created because reference genome is NULL")
  }

  if (!is.character(conf[['fillMissAllele']]) & !is.null(conf[['fillMissAllele']])) {
    stop("fillMissAllele option takes merged bam file as an argument!!!")
  } 				
  
  
  #==========================
  ### Create directiories ###
  #==========================
  ## Creating directories for data export
  if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
  }
  
  phased.store <- file.path(outputfolder, 'Phased')
  if (!file.exists(phased.store)) {
    dir.create(phased.store)
  }
  
  data.store <- file.path(outputfolder, 'data')
  if (!file.exists(data.store)) {
    dir.create(data.store)
  }
  
  browser.store <- file.path(outputfolder, 'browserFiles')
  if (!file.exists(browser.store)) {
    dir.create(browser.store)
  }
  
  vcf.store <- file.path(outputfolder, 'VCFfiles')
  if (!file.exists(vcf.store)) {
    dir.create(vcf.store)
  }
  
  singlecell.store <- file.path(outputfolder, 'SingleCellHaps')
  if (!file.exists(singlecell.store)) {
    dir.create(singlecell.store)
  }

  ## export haplotypes
  destination <- file.path(phased.store, 'phased_haps.txt')
  phased_haps.header <- matrix(c("sample", "cell", "chrom", "start", "end", "class", "hap1.cis.simil", "hap1.trans.simil", "hap2.cis.simil", "hap2.trans.simil"), nrow = 1)	
  write.table(data.frame(phased_haps.header), file=destination, row.names = F, quote = F, col.names = F, sep = "\t")	
  
  ## Make a copy of the conf file
  writeConfig(conf, configfile=file.path(outputfolder, 'StrandPhaseR.config'))
  
  ## Load BSgenome
  if (class(conf[['bsGenome']])!='BSgenome') {
    if (is.character(conf[['bsGenome']])) {
      suppressPackageStartupMessages(library(conf[['bsGenome']], character.only=T))
      conf[['bsGenome']] <- as.object(conf[['bsGenome']]) # replacing string by object
    }
  }
  
  ## Loading in list of SNV positions and locations of WC regions
  #snvs <- read.table(conf[['positions']], header=F)
  #snvs <- GRanges(seqnames=snvs$V1, IRanges(start=snvs$V2, end=snvs$V2))
  snvs <- vcf2ranges(vcfFile=conf[['positions']], genotypeField=1, chromosome=conf[['chromosomes']])
  #vcf <- read.vcfR(file = conf[['positions']], limit = 10000000, convertNA = TRUE, verbose = FALSE)
  #gt <- extract.gt(vcf)
  #hets <- is_het(gt)
  
  WC.regions <- read.table(conf[['WCregions']], header=F)
  #WC.regions <- read.table(conf[['WCregions']], header=F, sep = ":")
  WC.regions <- GRanges(seqnames=WC.regions$V1, IRanges(start=WC.regions$V2, end=WC.regions$V3), filename=as.character(WC.regions$V4))
  
  ## Parallelization ##
  if (conf[['numCPU']] > 1) {
    message("Using ",conf[['numCPU']]," CPUs")
    cl <- makeCluster(conf[['numCPU']])
    doParallel::registerDoParallel(cl)
  
    message("Phasing chromosomes ...", appendLF=F); ptm <- proc.time()
    temp <- foreach (chr = conf[['chromosomes']], .packages=c('StrandPhaseR')) %dopar% {
      chr <- as.character(chr)
    
      #Select chromosome of interest from the list
      snvs.chr <- snvs[seqnames(snvs) == chr]
      WCregions.chr <- WC.regions[seqnames(WC.regions) == chr]
      
      if (length(WCregions.chr) > 0 & length(snvs.chr) > 0) {
        seqlevels(snvs.chr) <- chr
        seqlevels(WCregions.chr) <- chr
        #snvs.chr <- keepSeqlevels(snvs.chr, chr, pruning.mode="coarse")
        #WCregions.chr <- keepSeqlevels(WCregions.chr, chr, pruning.mode="coarse")
        #snvs.chr <- keepSeqlevels(snvs.chr, chr)
        #WCregions.chr <- keepSeqlevels(WCregions.chr, chr)	
        
        tC <- tryCatch({
          phaseChromosome(inputfolder=inputfolder, outputfolder=outputfolder, positions=snvs.chr, WCregions=WCregions.chr, chromosome=chr, pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']], min.baseq=conf[['min.baseq']], num.iterations=conf[['num.iterations']], translateBases=conf[['translateBases']], fillMissAllele=conf[['fillMissAllele']], splitPhasedReads=conf[['splitPhasedReads']], compareSingleCells=conf[['compareSingleCells']], callBreaks=conf[['callBreaks']], exportVCF=conf[['exportVCF']], bsGenome=conf[['bsGenome']]) 
        }, error = function(err) {
          stop(chr,'\n',err)
        })
		
      } else {
        message("No SNVs or WC regions for a give chromosome ", chr)
      }		
    }
    stopCluster(cl)
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")  
    
  } else {
    
    for (chr in conf[['chromosomes']]) {
      chr <- as.character(chr)
      #Select chromosome of interest from the list
      snvs.chr <- snvs[seqnames(snvs) == chr]
      WCregions.chr <- WC.regions[seqnames(WC.regions) == chr]
      
      if (length(WCregions.chr) > 0 & length(snvs.chr) > 0) {
        seqlevels(snvs.chr) <- chr
        seqlevels(WCregions.chr) <- chr
        #snvs.chr <- keepSeqlevels(snvs.chr, chr, pruning.mode="coarse")
        #WCregions.chr <- keepSeqlevels(WCregions.chr, chr, pruning.mode="coarse")
        #snvs.chr <- keepSeqlevels(snvs.chr, chr)
        #WCregions.chr <- keepSeqlevels(WCregions.chr, chr)			
        
        phaseChromosome(inputfolder=inputfolder, outputfolder=outputfolder, positions=snvs.chr, WCregions=WCregions.chr, chromosome=chr, pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']], min.baseq=conf[['min.baseq']], num.iterations=conf[['num.iterations']], translateBases=conf[['translateBases']], fillMissAllele=conf[['fillMissAllele']], splitPhasedReads=conf[['splitPhasedReads']],  compareSingleCells=conf[['compareSingleCells']], callBreaks=conf[['callBreaks']], exportVCF=conf[['exportVCF']], bsGenome=conf[['bsGenome']]) 
	      
      }	else {
        message("No SNVs or WC regions for a give chromosome ", chr)	
      }
    }
  } 
} #end of wrapper function
