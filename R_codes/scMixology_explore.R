setwd('~/scLMM/sc_mixology/')
library(scran)
library(scater)
library(ggplot2)
library(DT)
library(Rtsne)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
#remotes::install_github(repo ='mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)
set.seed(5252)

library(Seurat)
library(SeuratData)
library(SeuratDisk)
#load("data/mRNAmix_qc.RData") # change to your own path
#load("data/9cellmix_qc.RData") # change to your own path

###########################################################################
############### Merging the 3cl data to be imported to python #############
###########################################################################

load("data/sincell_with_class.RData") ## 3 cell line data
#sce10x_qc contains the read counts after quality control processing from the 10x platform. 
#sce4_qc contains the read counts after quality control processing from the CEL-seq2 platform. 
#scedrop_qc_qc contains the read counts after quality control proessing from the Drop-seq platform.
data_list_3cl = list(sc_10x=sce_sc_10x_qc, 
                     CELseq2=sce_sc_CELseq2_qc, 
                     Dropseq=sce_sc_Dropseq_qc)

data_list_3cl = lapply(data_list_3cl, function(x) as.Seurat(x, counts = "counts", data = "logcounts"))
names(data_list_3cl) = c('sc_10X', 'CELseq2', 'Dropseq')
data_list_3cl = sapply(1:length(data_list_3cl), 
                       function(i) {data_list_3cl[[i]]$sample=names(data_list_3cl)[i]; data_list_3cl[[i]]}, simplify = F)
names(data_list_3cl) = c('sc_10X', 'CELseq2', 'Dropseq')

scMix_3cl_merged <- merge(data_list_3cl[[1]], c(data_list_3cl[[2]], data_list_3cl[[3]]),
                            add.cell.ids = names(data_list_3cl), 
                            project = "scMix_3cl", 
                            merge.data = TRUE)

table(scMix_3cl_merged$sample)
table(scMix_3cl_merged$cell_line)
table(scMix_3cl_merged$cell_line_demuxlet)
head(scMix_3cl_merged)



GetAssayData(scMix_3cl_merged, 'counts')[1:5, 1:5]
SaveH5Seurat(scMix_3cl_merged, filename = "~/scLMM/sc_mixology/scMix_3cl_merged.h5Seurat")
Convert("~/scLMM/sc_mixology/scMix_3cl_merged.h5Seurat", dest = "h5ad")

###########################################################################
############### Merging the 5cl data to be imported to python #############
###########################################################################
rm(list=ls())
load("data/sincell_with_class_5cl.RData") ## 5 cell line data
data_list_5cl = list(sc_10x=sce_sc_10x_5cl_qc, 
                     CELseq2_p1=sc_Celseq2_5cl_p1, 
                     CELseq2_p2=sc_Celseq2_5cl_p2,
                     CELseq2_p3=sc_Celseq2_5cl_p3)



data_list_5cl = lapply(data_list_5cl, function(x) logNormCounts(x))
data_list_5cl = lapply(data_list_5cl, function(x) as.Seurat(x, counts = "counts", data = "logcounts"))
data_list_5cl_names = names(data_list_5cl)
data_list_5cl = sapply(1:length(data_list_5cl), 
                       function(i) {data_list_5cl[[i]]$sample=names(data_list_5cl)[i]; data_list_5cl[[i]]}, simplify = F)
names(data_list_5cl) = data_list_5cl_names

scMix_5cl_merged <- merge(data_list_5cl[[1]], c(data_list_5cl[[2]], data_list_5cl[[3]], data_list_5cl[[4]]),
                          add.cell.ids = names(data_list_5cl), 
                          project = "scMix_5cl", 
                          merge.data = TRUE)



colnames(scMix_5cl_merged@meta.data)
sapply(1:ncol(scMix_5cl_merged@meta.data), function(i) table(scMix_5cl_merged@meta.data[[i]]))
table(scMix_5cl_merged$demuxlet_cls)
table(scMix_5cl_merged$cell_line)
table(scMix_5cl_merged$cell_line_demuxlet)
table(scMix_5cl_merged$batch)

head(scMix_5cl_merged)
table(scMix_5cl_merged$cell_line_demuxlet,scMix_5cl_merged$cell_line) 
table(scMix_5cl_merged$batch,scMix_5cl_merged$cell_line_demuxlet) 
table(scMix_5cl_merged$cell_line,scMix_5cl_merged$is_cell_control) 
table(scMix_5cl_merged$demuxlet_cls, scMix_5cl_merged$cell_line)


sum(is.na(scMix_5cl_merged$cell_line_demuxlet))
sum(is.na(scMix_5cl_merged$cell_line))
table(scMix_5cl_merged$cell_line_demuxlet[is.na(scMix_5cl_merged$cell_line)])
table(scMix_5cl_merged$batch[is.na(scMix_5cl_merged$cell_line)])
table(scMix_5cl_merged$demuxlet_cls[is.na(scMix_5cl_merged$cell_line)])
table(scMix_5cl_merged$[is.na(scMix_5cl_merged$cell_line)])

table(scMix_5cl_merged$cell_line_demuxlet, scMix_5cl_merged$cell_line)

sum(is.na(scMix_5cl_merged$cell_line_demuxlet))
sum(is.na(scMix_5cl_merged$cell_line))

GetAssayData(scMix_5cl_merged, 'counts')[1:5, 1:5]
SaveH5Seurat(scMix_5cl_merged, filename = "~/scLMM/sc_mixology/scMix_5cl_merged.h5Seurat")
Convert("~/scLMM/sc_mixology/scMix_5cl_merged.h5Seurat", dest = "h5ad")



head(scMix_5cl_merged)
lapply(data_list_5cl, function(x) table(x$cell_line))
##### Ground truth label information
#single cells trule label in colData: cell_line_demuxlet. 
#single cell mixtures the ground truth:  combination of three cell lines, column H1975, H2228 and HCC827. 
#paste(sce_SC1_qc$H1975,sce_SC1_qc$H2228,sce_SC1_qc$HCC827,sep="_")

#ground truth in RNA mixture is the proportion of RNA from each cell line, 
#column H2228_prop, H1975_prop and HCC827_prop
#paste(sce2_qc$H2228_prop,sce2_qc$H1975_prop,sce2_qc$HCC827_prop,sep="_")

table(sce_sc_10x_qc$cell_line)

sce_sc_10x_qc <- runPCA(sce_sc_10x_qc)
sce_sc_10x_qc.seur <- as.Seurat(sce_sc_10x_qc, counts = "counts", data = "logcounts")
head(sce_sc_10x_qc.seur@meta.data)
table(sce_sc_10x_qc.seur@meta.data$demuxlet_cls, sce_sc_10x_qc.seur@meta.data$cell_line)
table(sce_sc_10x_qc.seur@meta.data$cell_line, sce_sc_10x_qc.seur@meta.data$cell_line_demuxlet)
# gives the same results; but omits defaults provided in the last line
manno.seurat <- as.Seurat(manno)
Idents(manno.seurat) <- "cell_type1"
p1 <- DimPlot(manno.seurat, reduction = "PCA", group.by = "Source") + NoLegend()
p2 <- RidgePlot(manno.seurat, features = "ACTB", group.by = "Source")





