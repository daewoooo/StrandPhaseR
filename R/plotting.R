#' This function will take phased reads for each chromosome and plot them 
#' 
#' @param datapath containing split watson and crick reads per haplotype
#' @param plotChromosomes
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
		if (length(index)>0) {
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
		ggplt <- ggplot(long.dfplt) + geom_bar(aes(x=variable, y=value), stat="identity") + facet_wrap(~Filename, ncol=col.num)
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

plotSingleCellHaps <- function(data, file=NULL) {
  plt.df <- reshape2::melt(data, measure.vars = c('cons1.simil','cons2.simil'))  
  plt <- ggplot(plt.df, aes(y=value,x=start, color=variable)) + geom_step()  + facet_grid(CellID ~ .) + theme_bw() + theme(strip.text.y = element_text(angle=0), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + scale_color_manual(values = c("darkgoldenrod1", "dodgerblue2"))
  suppressMessages( ggsave(file, plot=plt, device="pdf", width=10, height=length(unique(plt.df$CellID))*0.5, limitsize=F) )
}



    
