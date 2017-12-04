#' Read StrandPhaseR configuration file
#'
#' Read an StrandPhaseR configuration file into a list structure. The configuration file has to be specified in INI format. R expressions can be used and will be evaluated.
#'
#' @param configfile Path to the configuration file
#' @author Aaron Taudt, David Porubsky
#' @export

readConfig <- function(configfile) {

	connection <- file(configfile) 
  Lines  <- readLines(connection) 
  close(connection) 

  Lines <- chartr("[]", "==", Lines) # change section headers 
	Lines <- gsub(" ", "", Lines) # no spaces

  connection <- textConnection(Lines) 
  data <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE, quote="") 
  close(connection) 
	names(data) <- c('argument','value','section')

  L <- data$argument == "" # location of section breaks 
  data$section <- data$value[which(L)[cumsum(L)]]
  data <- data[data$argument!="",]

  configlist <- list() 
	ToParse <- paste0("configlist$", data$argument, " <- ", data$value)
#   ToParse  <- paste0("configlist$", data$section, "$",  data$argument, " <- ", data$value) # with sections

  eval(parse(text=ToParse)) 

  return(configlist) 
} 

#' Write StrandPhaseR configuration file
#'
#' Write an StrandPhaseR configuration file from a list structure.
#'
#' @param conf A list structure with parameter values. Each entry will be written in one line.
#' @param configfile Filename of the outputfile.
#' @author Aaron Taudt, David Porubsky
#' @export

writeConfig <- function(conf, configfile) {

	## Printing function
	formatstring <- function(string) {
		if (is.character(string) & length(string)>1) {
			string <- paste0("c('",paste0(string,collapse="','"),"')")
		} else if (is.character(string) & length(string)==1) {
			string <- paste0("'",string,"'")
		} else if (is.numeric(string) & length(string)>1) {
			string <- paste0("c(",paste0(string,collapse=','),")")
		} else if (is.numeric(string) & length(string)==1) {
			string <- string
		} else if (is.null(string)) {
			string <- "NULL"
		}
		return(string)
	}
		
	f <- file(configfile, open='w')
	cat("#============== StrandPhaseR configuration file ===============#\n", file=f)
	cat("\n[General]\n", file=f)
	for (i1 in c('numCPU','chromosomes', 'pairedEndReads', 'min.mapq')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[StrandPhaseR]\n", file=f)
	for (i1 in c('positions', 'WCregions', 'min.baseq', 'num.iterations', 'translateBases', 'fillMissAllele', 'splitPhasedReads', 'compareSingleCells', 'callBreaks', 'exportVCF', 'bsGenome')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	close(f, type='w')
}
