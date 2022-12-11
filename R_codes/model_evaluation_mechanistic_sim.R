library(plyr)
library(stats)
library(ggpubr)
library(lme4) # load library
library(arm) # convenience functions for regression in R
library(dplyr)
library("pryr")
library(logr)
library(mice)
library(Seurat)
library(stringr)
library(MASS)
library(reshape2)
library(mltools)
library(data.table)
library('dbarts')
library(splatter)
library(scater)
library(purrr) # v. 0.3.4
library(broom) # v. 0.5.6
library(dplyr) # v. 1.0.0
library(ggplot2) # v. 3.3.1
library(Seurat)
library(SingleCellExperiment)


get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  #scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  scores          <- initial_data %*% invLoadings ## this second scaling is not necessary
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}


MODEL = 'LOG_LM'#'POIS_GLM' # 'NB_GLM', 'LOG_LM'
get_head<- function(df) df[1:4, 1:4]

# Empty for the sims means no grn
# 15 says how many time steps ~2^15cells x2 for this set (15-rp)
# The part of the name that you need is after the word empty—(gene number)



num_genes = 50
sampledata_files = list.files('simulated_data/cellmatrix_u2/', pattern = '*.txt', full.names = T, include.dirs = T)
sampledata_files_names = list.files('simulated_data/cellmatrix_u2/', pattern = '*.txt')
input_file = sampledata_files[grepl(pattern = paste0('empty--',num_genes,'-'), sampledata_files)]

file_number = 2
input_file = input_file[file_number]


#sampledata_list = lapply(sampledata_files, read.table)
data = read.table(input_file)
rownames(data) = paste0('gene', 1:nrow(data))
colnames(data) = paste0('cell', 1:ncol(data))

get_head(data)
dim(data)

seur = CreateSeuratObject(data)
dim(seur)

get_num_var_features <- function(seur){
  total_features = nrow(seur)
  if(total_features<10) return(1)
  else if (total_features < 20) return(2)
  else return(round(total_features)/10) 
} 

seur = FindVariableFeatures(seur)
# 2^15 cells - different number of genes
seur = ScaleData(seur)
seur <- SCTransform(seur,conserve.memory=F,verbose=T,
                    #variable.features.n = round(nrow(seur[['RNA']])/10), 
                    variable.features.n = get_num_var_features(seur), 
                    ## according to the paper scaling is not recommended prior to PCA:
                    ## https://www.biorxiv.org/content/10.1101/576827v2
                    do.scale = FALSE, ### default value 
                    do.center = TRUE) ### default value 

seur <- RunPCA(seur,verbose=T)
plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
num_PCs = 4#20
seur <- RunUMAP(seur,dims=1:num_PCs, reduction="pca")
df_umap <- data.frame(UMAP_1=getEmb(seur, 'umap')[,1], 
                      UMAP_2=getEmb(seur, 'umap')[,2],
                      library_size= seur$nCount_RNA, 
                      n_expressed=seur$nFeature_RNA )

dim(df_umap)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+
  geom_point(size=1)+theme_classic()+scale_color_viridis(direction = -1)+
  ggtitle(paste('gene number:', num_genes, ' #PCs:', num_PCs, '\n', 
                substr(sampledata_files_names[file_number], start=1,stop = 60), '\nfile-number: ', file_number))

seur <- RunTSNE(seur,dims=1:num_PCs, reduction="pca")
df_tsne <- data.frame(tSNE_1=getEmb(seur, 'tsne')[,1], 
                      tSNE_2=getEmb(seur, 'tsne')[,2],
                      library_size= seur$nCount_RNA, 
                      n_expressed=seur$nFeature_RNA )

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+
  geom_point(size=1)+theme_classic()+scale_color_viridis(direction = -1)+
  ggtitle(paste('gene number:', num_genes, ' #PCs:', num_PCs, '\n', 
                substr(sampledata_files_names[file_number], start=1,stop = 70)))

# all genes are on and tehy are all accesible - the model is not biased at all

rm(input_file)
