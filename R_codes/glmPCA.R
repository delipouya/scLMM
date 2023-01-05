library(glmpca)
library(splatter)
library(scater)
library(purrr) # v. 0.3.4
library(broom) # v. 0.5.6
library(dplyr) # v. 1.0.0
library(ggplot2) # v. 3.3.1
library(Seurat)
library(SingleCellExperiment)
get_head<- function(df) df[1:4, 1:4]

########################################################################
############ Simple dataset offered by the glmPCA package ############
#create a simple dataset with two clusters
mu<-rep(c(.5,3),each=10)
mu<-matrix(exp(rnorm(100*20)),nrow=100)
mu[,1:10]<-mu[,1:10]*exp(rnorm(100))
clust<-rep(c("red","black"),each=10)
Y<-matrix(rpois(prod(dim(mu)),mu),nrow=nrow(mu))
dim(Y)
############################################################

sim = readRDS('simulated_data/simulated_data_3groups_2batch.rds')

counts(sim)[1:5, 1:5] ## count data is simulated
head(rowData(sim)) ### can be used to compare with coefficients and the residual based factors
head(colData(sim)) 
names(assays(sim))
assays(sim)$CellMeans[1:5, 1:5]

plotUMAP(sim, shape_by = "Batch", colour_by = "Group")
plotUMAP(sim, shape_by = "Group", colour_by = "Batch")
plotPCA(sim, shape_by = "Batch", colour_by = "Group")
plotPCA(sim, shape_by = "Group", colour_by = "Batch")

gene_info = as.data.frame(rowData(sim))
meta_data = as.data.frame(colData(sim))
head(meta_data)
head(gene_info)
dim(sim)
table(meta_data$Batch)
table(meta_data$Group)

Y = counts(sim)
get_head(Y)
gene_idx = rowSums(counts(sim))!=0

X = as.matrix(ifelse(meta_data$Group=='Group1', 0, 
                     ifelse(meta_data$Group=='Group2', 1, 2)))
X = as.matrix(ifelse(meta_data$Batch=='Batch1', 0, 1))

Y = Y[gene_idx,]
gene_info = gene_info[gene_idx,]
dim(Y)
dim(gene_info)

r = 5
#visualize the latent structure
res<-glmpca(Y=Y, X=X, L=r, verbose=T)
factors<-res$factors
plot(factors[,1],factors[,2],col=clust,pch=19)
head(res$coefX)


score_df = data.frame(cbind(res$factors, meta_data))
loading_df = as.data.frame(cbind(res$loadings, gene_info))

coefficient_df = as.data.frame(res$coefX)
head(coefficient_df)
colnames(coefficient_df) = c('intercepts', 'cov1')


factor_names =  paste0('factor_', 1:r)
colnames(score_df)[1:r] = factor_names
colnames(loading_df)[1:r] = factor_names

ggplot(score_df, aes(x=factor_1, y=factor_2, color=Batch))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_3, color=Batch))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_4, color=Batch))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_5, color=Batch))+geom_point()+theme_classic()

ggplot(score_df, aes(x=factor_1, y=factor_2, color=Group))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_3, color=Group))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_4, color=Group))+geom_point()+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_5, color=Group))+geom_point()+theme_classic()

ggplot(loading_df, aes(x=factor_1, y=DEFacGroup1, color=factor_2))+geom_point()+theme_classic()
ggplot(loading_df, aes(x=factor_2, y=BatchFacBatch2, color=GeneMean))+geom_point()+theme_classic()
ggplot(loading_df, aes(x=factor_2, y=BatchFacBatch1, color=GeneMean))+geom_point()+theme_classic()
head(loading_df)

cor_mat = cor(loading_df[,c(1:5,7,9:ncol(loading_df))])
cor_mat = cor(V[,c(1:3,5,7:ncol(V))])
pheatmap(cor_mat[rownames(cor_mat) %in% factor_names, !rownames(cor_mat) %in% factor_names])
pheatmap(res_sim$Sf)
#############################################################

###################################################################
########## applying the glmPCA to the rat samples ############## 
##raw counts - importing the data and quality control

library(Seurat)
library(SeuratDisk)
library(SeuratData)
rat_h5seur = LoadH5Seurat("~/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5seurat")


dat <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
sample_info = sapply(strsplit(colnames(dat), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
strain_info = sapply(strsplit(colnames(dat), '_'), '[[', 1)
table(sample_info)

meta_data = data.frame(cell_id = colnames(dat), sample=sample_info, 
                       strain=strain_info, cluster=rat_h5seur$cluster)

dim(dat) #[1] 32883 23036
get_head(dat)
ncells <- rowSums(dat > 0)
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))

genes_2keep = ncells >= 20
Y <- dat[genes_2keep, ]
genelist = rownames(Y)
rm(dat)
gc()
get_head(Y)


#Y = as.matrix(ifelse(meta_data$Group=='Group1', 0, ifelse(meta_data$Group=='Group2', 1, 2)))
X = as.matrix(ifelse(strain_info=='DA', 0, 1))
table(X)
###### Supervided SVD tutorial
r=10 ## number of factors
dim(X)
dim(Y)
head(Y)
colnames(X) = paste0('cov_', 1:ncol(X)) # strain

start_time <- Sys.time()
res<-glmpca(Y=Y, X=X, L=r, verbose=T)
end_time <- Sys.time()
print(paste0('run time: ', end_time - start_time))
res <- readRDS('~/scLMM/glmPCA_ratTLH_strain_cov_result.rds')
#SupPCA_res = readRDS('~/scLMM/glmPCA_ratTLH_strain_cov_result.rds')


score_df = data.frame(cbind(res$factors, meta_data))
loading_df = as.data.frame(cbind(res$loadings))

coefficient_df = as.data.frame(res$coefX)
head(coefficient_df)
colnames(coefficient_df) = c('intercepts', 'cov1')


factor_names =  paste0('factor_', 1:r)
colnames(score_df)[1:r] = factor_names
colnames(loading_df)[1:r] = factor_names

ggplot(score_df, aes(x=factor_1, y=factor_2, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_2, color=sample))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_2, color=strain))+geom_point(size=1)+theme_classic()

ggplot(score_df, aes(x=factor_1, y=factor_3, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_3, color=sample))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_3, color=strain))+geom_point(size=1)+theme_classic()

ggplot(score_df, aes(x=factor_1, y=factor_4, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_4, color=sample))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_4, color=strain))+geom_point(size=1)+theme_classic()

ggplot(score_df, aes(x=factor_1, y=factor_5, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_5, color=sample))+geom_point(size=1)+theme_classic()
ggplot(score_df, aes(x=factor_1, y=factor_5, color=strain))+geom_point(size=1)+theme_classic()


ggplot(score_df, aes(x=factor_1, y=factor_10, color=cluster))+geom_point()+theme_classic()
head(loading_df)


