
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
pkg = "xopen"
pkg = 'devtools'
pkg = c('textshaping', 'ragg', 'pkgdown')
pkg = 'dbplyr'
pkg = 'combinat'
pkg = 'MCMCprecision'#'fontquiver'

i = 10
pkgs = c('numDeriv', 'ggplotify','interactiveDisplayBase', 'AnnotationHub', 
         'annotate', 'GSEABase', 'genefilter', 'ExperimentHub', 'scran')
pkg = pkgs[i]
pkg = 'ggtree'
pkg = 'singleCellTK'
BiocManager::install(pkg, version = '3.15', lib =.libPaths()[1], force = TRUE)

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

BiocManager::install(pkg, version = '3.15', lib =.libPaths()[1], force = TRUE)
remove.packages(c(pkg), lib = .libPaths()[1])
remove.packages(c(pkg), lib = .libPaths()[2])
remove.packages(c(pkg), lib = .libPaths()[3])
remove.packages(c(pkg), lib = .libPaths()[4])


#https://github.com/astamm/nloptr/issues/40
#https://stackoverflow.com/questions/37425509/trouble-installing-nloptr-package-on-r-3-3-0

#sudo apt-get install libnlopt-dev

library(tools)
db <- CRAN_package_db()
revdeps <- package_dependencies("Seurat", db=db, recursive=TRUE, reverse=TRUE)



install.packages("tensorflow")
base64enc
tfruns
