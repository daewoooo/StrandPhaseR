#' This function will take phased reads for each chromosome and plot them 
#' 
#' @param datapath containing split watson and crick reads per haplotype
#' @param plotChromosomes ...
#' @param file name of the file to save produced plots
#' 
#' @author David Porubsky
#'
plotPhasedReads <- function(datapath, plotChromosomes = NULL, file=NULL) {
  
  RData.objects <- list.files(datapath, pattern="_reads.RData$", full=T)
  
  data <- get(load(RData.objects))  
  
  hap1 <- data$hap1
  hap2 <- data$hap2
  
  #get coverage per genomic position
  rle.hap1.cov <- coverage(hap1)
  hap1.cov <- runValue(rle.hap1.cov) # get values of coverage per genomic site 
  hap1.cov <- unlist(hap1.cov, use.names = FALSE) # transfer coverage values into vector
  #hap1.cov <- hap1.cov + 1 # add value 1 to coverage to reach join line across whole chromosome => baseline coverage is one	
  hap1.cov.ranges <- ranges(rle.hap1.cov) # get genomic ranges with certain value of coverage 
  hap1.cov.ranges <- unlist(hap1.cov.ranges, use.names = FALSE) 
  df.hap1Cov <- data.frame(start=start(hap1.cov.ranges), end=end(hap1.cov.ranges), cov=hap1.cov)
  
  rle.hap2.cov <- coverage(hap2)
  hap2.cov <- runValue(rle.hap2.cov)
  hap2.cov <- unlist(hap2.cov, use.names = FALSE)
  #hap2.cov <- hap2.cov + 1	
  hap2.cov.ranges <- ranges(rle.hap2.cov)
  hap2.cov.ranges <- unlist(hap2.cov.ranges, use.names = FALSE)
  df.hap2Cov <- data.frame(start=start(hap2.cov.ranges), end=end(hap2.cov.ranges), cov=hap2.cov)
  
  ggplt <- ggplot(df.hap1Cov, aes(x=start, y=cov))
  ggplt <- ggplt + geom_step(col='red')
  ggplt <- ggplt + geom_step(data=df.hap2Cov, aes(x=start, y=-cov), col = 'blue', inherit.aes = F)
  
  return(ggplt)
}


#' Load phased data from RData files
#'
#' @param files A list of files that contain phased data.
#' @return A list() containing all loaded phased data.
#' @author David Porubsky
#'
loadGRangesFromFiles <- function(files) {

	file.list <- files
	if (is.character(file.list)) {
		message("Loading phased data from files ...", appendLF=F); ptm <- proc.time()
		phasedlist <- list()
		for (phasedChrom in file.list) {
			tC <- tryCatch({
				phasedlist[[phasedChrom]] <- get(load(phasedChrom))
			}, error = function(err) {
				stop(phasedChrom,'\n',err)
			})
			#if (class(mlist[[modelfile]])!='GRanges') {
			#	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
			#	stop("File ",modelfile," does not contain a GRanges object.")
			#}
		}
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
		return(phasedlist)
	} else if (is.list(file.list)) {
		index <- which(unlist(lapply(gr.list, function(hmm) { class(hmm)!='GRanges' })))
		if (length(index) > 0) {
			stop("The following list entries do not contain GRanges objects: ", paste(index, collapse=' '))
		}
		return(gr.list)
	}
}


#' This function will take phased data for each chromosome and plot density of phased SNVs 
#' 
#' @param datapath containing split watson and crick reads per haplotype
#' @param perChromosome if TRUE function will plot density for every chromosome, otherwise one summary density plot is exported
#' @param file name of the file to save produced plots (add appropriate file extension)
#' 
#' @author David Porubsky
#'
plotHapDensity <- function(datapath, perChromosome=FALSE, file=NULL, HighQual=FALSE) {

	phased.data <- list.files(datapath, pattern="_phased_hap", full=T)	

	files.dist <- list()
	for (i in 1:length(phased.data)) {

		phasedHap <- read.table(phased.data[i], header=T)
		if (HighQual) {
			qual.factor <- (1-phasedHap$ent)*phasedHap$cov
			pos <- phasedHap$pos[qual.factor>=2]
		} else {
			pos <- phasedHap$pos
		}

		dist <- sort(diff(pos)-1)
		categ <- c(1001,1501,2001,2501,3001,3501,4001,4501,5001)
		categ.intervals <- findInterval(dist, categ)
		categ.counts <- table(categ.intervals)
		files.dist[[basename(phased.data[i])]] <- as.vector(categ.counts)
	}

	categories <- c("0-1000", "1001-1500", "1501-2000", "2001-2500", "2501-3000", "3001-3500", "3501-4000", "4001-4500", "4501-5000", ">5000")
	phased.all.chr <- do.call(rbind, files.dist)

	if (perChromosome) {
		dfplt <- as.data.frame(phased.all.chr)
		names(dfplt) <- categories
		dfplt$Filename <- rownames(dfplt)
		col.num <- floor(sqrt(length(rownames(dfplt))))
		long.dfplt <- melt(dfplt)
		ggplt <- ggplot(long.dfplt) + 
		  geom_bar(aes(x=variable, y=value), stat="identity") + 
		  facet_wrap(~Filename, ncol=col.num)
		if (!is.null(file)) {
			filetype <- strsplit(file, "\\.")[[1]][2]
			ggsave(file, ggplt, width=col.num*2 ,height=col.num*2, limitsize=FALSE, device=filetype)
		}
	} else {
		dist.sum <- colSums(phased.all.chr)
		dfplt <- data.frame(Categ=categories, value=dist.sum)
		ggplt <- ggplot(dfplt) + geom_bar(aes(x=Categ, y=value), stat="identity")
		if (!is.null(file)) {
			filetype <- strsplit(file, "\\.")[[1]][2]
			ggsave(file, ggplt, width=length(categories) ,height=length(categories)/2, device=filetype)
		}
	}
	if(is.null(file)) { return(ggplt) }
}		
		
#' This function will take phased data for each chromosome and plot density of phased SNVs 
#' 
#' @param data A \code{list} output from the function \code{\link{compareSingleCellHaps}}
#' @param file name of the file to save produced plots (add appropriate file extension)
#' 
#' @author David Porubsky
#' @export
#'
plotSingleCellHaps <- function(data, file=NULL) {
  plt.df <- reshape2::melt(data, measure.vars = c('cons1.simil','cons2.simil'))  
  plt <- ggplot(plt.df, aes(y=value,x=start, color=variable)) + 
    geom_step()  + 
    facet_grid(CellID ~ .) + 
    theme_bw() + 
    theme(strip.text.y = element_text(angle=0), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
    scale_color_manual(values = c("darkgoldenrod1", "dodgerblue2"))
  suppressMessages( ggsave(file, plot=plt, device="pdf", width=10, height=length(unique(plt.df$CellID))*0.5, limitsize=F) )
}

#' Plot inherited haplotype segments as ideogram.
#' 
#' This function takes as na input a \code{\link{GRanges-class}} object with genomic positions of each
#' inherited haplotype block and project them as whole-genome ideogram.
#' 
#' @param hapSegm.gr A \code{\link{GRanges-class}} object with regions of each inherited haplotype block
#' @param layout Set to 'vertical' or 'horizonal' based on desired layout of an resultant ideogram.
#' @importFrom dplyr "%>%"
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#' 
plotHaploSegs <- function(hapSegm.gr=NULL, layout='vertical') {
  ## Remove ranges with NA values
  hapSegm.gr <- hapSegm.gr[!is.na(hapSegm.gr$match)]
  ## Prepare data for plotting
  plt.df <- as.data.frame(hapSegm.gr)
  ## Create unique ID
  plt.df$ID <- paste0(plt.df$match, ".", plt.df$inherited)
  ## Set levels for homolog plotting
  homolog.ids <- unique(plt.df$inherited)
  plt.df$ymin <- 0
  plt.df$ymax <- 0
  plt.df$ymin[plt.df$inherited == homolog.ids[1]] <- 0
  plt.df$ymax[plt.df$inherited == homolog.ids[1]] <- 2
  plt.df$ymin[plt.df$inherited == homolog.ids[2]] <- 3
  plt.df$ymax[plt.df$inherited == homolog.ids[2]] <- 5
  
  ## Set plotting themes
  theme_horizontal <- theme(legend.position ="top", 
                            axis.text.y=element_blank(), 
                            strip.text.y = element_text(angle = 180),
                            axis.ticks.y = element_blank(),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank())
  
  theme_vertical <- theme(legend.position ="top",
                          axis.line = element_blank(),
                          axis.text.x=element_blank(), 
                          axis.ticks.x=element_blank(),   
                          strip.text.y = element_text(angle = 180),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(plt.df$end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(plt.df$seqnames))
  
  ## Count number of breakpoints per homolog
  break.count <- plt.df %>% 
    group_by(inherited, seqnames) %>% 
    summarise(count=n() - 1) %>% 
    group_by(inherited) %>% 
    summarise(count=sum(count))
  ## Prepare breakpoint annotation
  break.count.annot <- paste0(break.count$inherited, ": ", break.count$count)
  #break.count.annot <- as.data.frame(break.count.annot)
  #break.count.annot$level <- 1:nrow(break.count.annot)
  break.count.annot <- paste0(break.count.annot, collapse = "\n")
  
  if (layout == 'horizontal') {
    plt <- ggplot(plt.df) + 
      geom_rect(aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=ID)) + 
      facet_grid(seqnames ~ ., switch = 'y') + 
      scale_fill_manual(values=c("cadetblue3","coral", "dodgerblue4", "firebrick"), name=break.count.annot) +
      scale_x_continuous(breaks = breaks, labels = labels) +
      theme_horizontal
  } else if (layout == 'vertical') {
    plt <- ggplot(plt.df) + 
      geom_rect(aes(xmin=ymin, xmax=ymax, ymin=start, ymax=end, fill=ID)) + 
      facet_grid(. ~ seqnames, switch = 'x') + 
      scale_fill_manual(values=c("cadetblue3","coral", "dodgerblue4", "firebrick"), name=break.count.annot) + 
      scale_y_continuous(breaks = breaks, labels = labels) +
      theme_vertical
  } else {
    message("Unsupported layout!!! Please choose from 'horizontal' or 'vertical'.")
  }
  return(plt)
}


    
