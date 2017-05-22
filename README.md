<img src="https://github.com/daewoooo/StrandPhaseR/raw/master/StrandPhaseR_logo.png" />
=========================================================================

# StrandPhaseR
R Package for phasing of single cell Strand-seq data 

Collaborators: David Porubsky, Ashley D. Sanders

#Installation

### Bioconductor version (not available yet)
Under the development.

### Development version from Github
To install the development version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>=3.2.0) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   	install.packages("devtools")
	source("http://bioconductor.org/biocLite.R")
	biocLite("GenomicRanges")
	biocLite("GenomicAlignments")
	library(devtools)
	install_github("daewoooo/StrandPhaseR")
	Or alternatively if the above line doesn't work:
	install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master")

### How to use BreakPointR

1. Start Rstudio
2. Load StrandPhaseR package: 	library('StrandPhaseR')
3. Run StrandPhaseR: 	strandPhaseR(inputfolder = <input_folder>, outputfolder = <output_folder>, positions = <SNV_positions>, WCregions = <haplotypeInformative_WCregions>, chromosomes = <chromosomes2analyze>) or use config file as follows strandPhaseR(inputfolder = <input_folder>, outputfolder = <output_folder>, positions = <SNV_positions>, WCregions = <haplotypeInformative_WCregions>, configfile = <config_file>)

### Report Errors

If you encounter errors of any kind, please report an [issue here](https://github.com/daewoooo/StrandPhaseR/issues/new).
