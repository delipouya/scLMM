
# create a list of all installed packages
ip <- as.data.frame(installed.packages())
head(ip)
# if you use MRO, make sure that no packages in this library will be removed
ip <- subset(ip, !grepl("MRO", ip$LibPath))
# we don't want to remove base or recommended packages either\
ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# determine the library where the packages are installed
path.lib <- unique(ip$LibPath)
# create a vector with all the names of the packages you want to remove
pkgs.to.remove <- ip[,1]
head(pkgs.to.remove)
# remove the packages
sapply(pkgs.to.remove, remove.packages, lib = path.lib)


devtools::install_github('RfastOfficial/Rfast')
library(BiocManager)
BiocManager::install(pkgs = ,version = '3.15',lib =.libPaths()[1], force = TRUE)
library(devtools)
library(nebula)
BiocManager::install()
library(RcppZiggurat)

packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_2.0.0.tar.gz"
packageurl<-'https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.0.4.tar.gz'
packageurl<-'https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.0.tar.gz'
packageurl<-'https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.5.tar.gz'
packageurl<-'https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.1.0.tar.gz'
packageurl<-'https://cran.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.3-2.tar.gz'
install.packages(packageurl, repos=NULL, type="source")
install.packages('Seurat')

.libPaths()
pkg = "SeuratDisk" #Seurat
pkg = 'roxygen2' 
pkg = c('textshaping', 'ragg', 'pkgdown')
pkg = 'lazyeval'
pkg = 'DelayedArray'
pkg = 'DelayedMatrixStats'#'fontquiver'

i = 10
pkgs = c('numDeriv', 'ggplotify','interactiveDisplayBase', 'AnnotationHub', 
         'annotate', 'GSEABase', 'genefilter', 'ExperimentHub', 'scran')
pkg = pkgs[i]
pkg = 'ggtree'
pkg = 'pscl'
install.packages(pkg)
pkg='graphlayouts'
#BiocManager::install(pkg, version = '3.15', lib =.libPaths()[1], force = TRUE)
BiocManager::install(pkg, version = '3.18', lib =.libPaths()[1], force = TRUE)
#testthat, xopen
library('pscl')
library(pkg)
BiocManager::install("ggtree")

install.packages("RcppZiggurat")
library(RcppZiggurat)
devtools::install_github("lc5415/HDATDS")
library(reticulate)
py_install("numpy")

devtools::install_github("lhe17/nebula")
install.packages("nebula", repos="http://R-Forge.R-project.org")

install.packages("Rfast")
library(Rfast)
BiocManager::install(pkg)

BiocManager::install(pkg, version = '3.18', lib =.libPaths()[1], force = TRUE)
remove.packages(c(pkg), lib = .libPaths()[1])
remove.packages(c(pkg), lib = .libPaths()[2])
remove.packages(c(pkg), lib = .libPaths()[3])
remove.packages(c(pkg), lib = .libPaths()[4])

Error: package ‘GenomicAlignments’ was installed before R 4.0.0: please re-install it
Error: package or namespace load failed for ‘dbplyr’:
Error: package ‘GSEABase’ was installed before R 4.0.0: please re-install it
Error: package ‘progress’ was installed before R 4.0.0: please re-install it
Error: package ‘AnnotationHub’ 2.18.0 was found, but >= 3.3.6 is required by ‘ExperimentHub’
Error: package ‘biomaRt’ was installed before R 4.0.0: please re-install it
Error: package ‘GenomicFeatures’ 1.38.1 was found, but >= 1.49.6 is required by ‘ensembldb’

package ‘DelayedArray’ was installed before R 4.0.0: please re-install it




#https://github.com/astamm/nloptr/issues/40
#https://stackoverflow.com/questions/37425509/trouble-installing-nloptr-package-on-r-3-3-0

#sudo apt-get install libnlopt-dev

library(tools)
db <- CRAN_package_db()
revdeps <- package_dependencies("Seurat", db=db, recursive=TRUE, reverse=TRUE)


# check your package library path 
.libPaths()

# grab old packages names
old_packages <- installed.packages(lib.loc =  "/usr/lib/R/library")

old_packages <- as.data.frame(old_packages)
list.of.packages <- unlist(old_packages$Package)

# remove old packages 
remove.packages( installed.packages( priority = "NA" )[,1] )

# reinstall all packages 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
