source('~/RatLiver/Codes/Functions.R')
source('~/RatLiver/Codes/FactorAnalysisUtils.R')
Initialize()

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

get_head<- function(df) df[1:4, 1:4]
library(SuperPCA)

##### evaluating the LM models on the simulated data using the splatter

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


X = t(counts(sim))
get_head(X)
Y = as.matrix(ifelse(meta_data$Group=='Group1', 0, ifelse(meta_data$Group=='Group2', 1, 2)))
Y = as.matrix(ifelse(meta_data$Batch=='Batch1', 0, 1))
table(Y)
###### Supervided SVD tutorial
r=3 ## number of factors
dim(X)
dim(Y)
head(Y)
colnames(Y) = paste0('cov_', 1:ncol(Y))
rownames(Y) = rownames(X)

Yc <- scale(Y,center=TRUE,scale=FALSE)
Xc <- scale(X,center=TRUE,scale=FALSE)

start_time <- Sys.time()
SupPCA_res = SupPCA(Yc,Xc,r)
#SupSFPCA_res = SupSFPCA(Yc,Xc,r) # error in glmnet x should be a matrix with 2 or more columns
end_time <- Sys.time()

res_sim = SupPCA_res
print(paste0('run time: ', end_time - start_time))
#saveRDS(res_sim, '~/scLMM/supsvd_sim_output.rds')
#res_sim = readRDS('~/scLMM/supsvd_sim_output.rds')
U = data.frame(cbind(res_sim$U, meta_data))
V = as.data.frame(cbind(res_sim$V, gene_info))
rownames(V) = colnames(X)
B = as.data.frame(res_sim$B)
print(B)
factor_names =  paste0('factor_', 1:r)
colnames(U)[1:r] = factor_names
colnames(V)[1:r] = factor_names
head(U)
head(V)
ggplot(U, aes(x=factor_1, y=factor_2, color=Batch))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_3, color=Batch))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_2, color=Group))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_3, color=Group))+geom_point()+theme_classic()

ggplot(V, aes(x=factor_1, y=DEFacGroup1, color=factor_2))+geom_point()+theme_classic()
ggplot(V, aes(x=factor_2, y=BatchFacBatch2, color=GeneMean))+geom_point()+theme_classic()

head(V)
cor_mat = cor(V[,c(1:3,5,7,8:ncol(V))])
cor_mat = cor(V[,c(1:3,5,7:ncol(V))])
pheatmap(cor_mat[rownames(cor_mat) %in% factor_names, !rownames(cor_mat) %in% factor_names])
pheatmap(res_sim$Sf)


####################################################################
########  Applying the Supervised SVD to the RAW rat liver data ########
####################################################################
##raw counts - importing the data and quality control
dat <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
sample_info = sapply(strsplit(colnames(dat), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
strain_info = sapply(strsplit(colnames(dat), '_'), '[[', 1)
table(sample_info)

meta_data = data.frame(cell_id = colnames(dat), sample=sample_info, strain=strain_info)
  
dim(dat) #[1] 32883 23036
get_head(dat)
ncells <- rowSums(dat > 0)
hist(ncells)
hist(log2(ncells))
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))

genes_2keep = ncells >= 20
X <- as.matrix(t(dat[genes_2keep, ]))
genelist = rownames(X)
rm(dat)
gc()
get_head(X)


#Y = as.matrix(ifelse(meta_data$Group=='Group1', 0, ifelse(meta_data$Group=='Group2', 1, 2)))
Y = as.matrix(ifelse(strain_info=='DA', 0, 1))
table(Y)
###### Supervided SVD tutorial
r=10 ## number of factors
dim(X)
dim(Y)
head(Y)
colnames(Y) = paste0('cov_', 1:ncol(Y))
rownames(Y) = rownames(X)

Yc <- scale(Y,center=TRUE,scale=FALSE)
Xc <- scale(X,center=TRUE,scale=FALSE)

start_time <- Sys.time()
#SupPCA_res = SupPCA(Yc,Xc,r)
#SupSFPCA_res = SupSFPCA(Yc,Xc,r) # error in glmnet x should be a matrix with 2 or more columns
end_time <- Sys.time()

print(paste0('run time: ', end_time - start_time))
SupPCA_res = readRDS('~/scLMM/supsvd_ratTLH_output.rds')
SupPCA_res$U

U = data.frame(cbind(SupPCA_res$U, meta_data))

V = as.data.frame(cbind(SupPCA_res$V))
rownames(V) = colnames(X)
B = as.data.frame(SupPCA_res$B)
print(B)
names(which.max(abs(B))) ## best match with the covariate of your interest

factor_names =  paste0('factor_', 1:r)
colnames(U)[1:r] = factor_names
colnames(V)[1:r] = factor_names
head(U)
head(V)

sapply(1:10, function(i) U[,i]<-as.numeric(U[,i]))
ggplot(U, aes(x=factor_1, y=factor_2, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_3, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_4, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_5, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_6, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_7, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_8, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_9, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_10, color=strain))+geom_point()+theme_classic()
head(V[order(V$factor_9, decreasing = T),], 20)

head(V[order(V$factor_2, decreasing = F),], 20)


head(U)
library(Seurat)
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
rm(your_scRNAseq_data_object); gc()
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'LEW-1', 'LEW-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
processed_df = data.frame(cluster=merged_samples$cluster, 
                          sample_name=merged_samples$sample_name, 
                          strain=merged_samples$strain, 
                          cell_id=colnames(merged_samples))

sample_info = sapply(strsplit(processed_df$cell_id, '_'), function(x) paste0(x[2], '_' ,x[3]))
pure_cell_id = sapply(strsplit(processed_df$cell_id, '_'), function(x) paste0(x[length(x)]))
sample_info[sample_info=='DA_M'] = 'DA_02'
sample_info[sample_info=='Lew_01'] = 'LEW_01'
sample_info[sample_info=='Lew_02'] = 'LEW_02'
table(sample_info)

##### channging the cell ID of the processed data
processed_df$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
head(processed_df)

######### changing the cell_iD of U
sample_info = sapply(strsplit(U$cell_id, '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
table(sample_info)
pure_cell_id = sapply(strsplit(U$cell_id, '_'), function(x) paste0(x[length(x)]))
U$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
plot(U[,1],U[,2])



umap_df = data.frame(merged_samples[["umap"]]@cell.embeddings)
umap_df$cell_id = processed_df$cell_id
umap_df$cell_id_2 = processed_df$cell_id_2

sum(!U$cell_id_2 %in% umap_df$cell_id_2)
sum(!umap_df$cell_id_2 %in% U$cell_id_2)
head(umap_df$cell_id_2[!umap_df$cell_id_2 %in% U$cell_id_2])
merged_df = merge(U, umap_df, by.x='cell_id_2', by.y='cell_id_2')
dim(U)
dim(merged_df)
head(merged_df)
Var_to_check = 'factor_1'
ggplot(merged_df, aes(UMAP_1, UMAP_2, color=factor_1))+
  geom_point()+theme_bw()




#######################################################################
############################################################################
########  Applying the Supervised SVD to the normalized rat liver data ########
############################################################################
#######################################################################
dat_norm = merged_samples
table(dat_norm$sample_name)
table(merged_samples$strain)

sample_info = sapply(strsplit(colnames(dat_norm), '_'), function(x) paste0(x[2], '_' ,x[3]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
sample_info[sample_info=='DA_M'] = 'DA_02'
sample_info[sample_info=='Lew_01'] = 'LEW_01'
sample_info[sample_info=='Lew_02'] = 'LEW_02'

strain_info = sapply(strsplit(colnames(dat_norm), '_'), '[[', 1)
table(sample_info)

meta_data = data.frame(cell_id = colnames(dat_norm), sample=sample_info, strain=strain_info, cluster=dat_norm$cluster)
head(meta_data)

X <- as.matrix(t(GetAssayData(dat_norm)))
genelist = rownames(X)
rm(dat_norm)
gc()
get_head(X)

library('plm')
c("LEW", "DA_01", "LEW_01")
#Y = as.matrix(ifelse(meta_data$Group=='Group1', 0, ifelse(meta_data$Group=='Group2', 1, 2)))
Y = as.matrix(data.frame(DA_01=ifelse(meta_data$sample=='DA_01', 0,1), 
                         LEW_01=ifelse(meta_data$sample=='LEW_01', 0,1),
                         DA_02=ifelse(meta_data$sample=='DA_02', 0,1) )
              ) #strain=ifelse(meta_data$strain=='DA', 0, 1), 
head(Y)
table(Y)

###### Supervided SVD tutorial
r=10 ## number of factors
dim(X)
dim(Y)
head(Y)
colnames(Y) = paste0('cov_', 1:ncol(Y))
rownames(Y) = rownames(X)

Yc <- scale(Y,center=TRUE,scale=FALSE)
head(Yc)
detect.lindep(Y)
detect.lindep(Yc)

Xc <- scale(X,center=TRUE,scale=FALSE)

start_time <- Sys.time()
SupPCA_res = SupPCA(Yc,Xc,r)
#SupSFPCA_res = SupSFPCA(Yc,Xc,r) # error in glmnet x should be a matrix with 2 or more columns
end_time <- Sys.time()
print(paste0('time is: ', end_time-start_time))
#saveRDS(SupPCA_res, '~/scLMM/supsvd_ratTLH_scaled_output.rds')
#saveRDS(SupPCA_res, '~/scLMM/supsvd_ratTLH_scaled_strain.cov.rds')

SupPCA_res <- readRDS('~/scLMM/supsvd_ratTLH_scaled_sample.cov.rds')

U = data.frame(cbind(SupPCA_res$U, meta_data))
V = as.data.frame(cbind(SupPCA_res$V))
rownames(V) = colnames(X)
B = as.data.frame(SupPCA_res$B)
print(B)
names(which.max(abs(B[1,]))) ## best match with the covariate of your interest
names(which.max(abs(B[2,]))) ## best match with the covariate of your interest
names(which.max(abs(B[3,]))) ## best match with the covariate of your interest

factor_names =  paste0('factor_', 1:r)
colnames(U)[1:r] = factor_names
colnames(V)[1:r] = factor_names
head(U)
head(V)

sapply(1:10, function(i) U[,i]<-as.numeric(U[,i]))
ggplot(U, aes(x=factor_1, y=factor_2, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_2, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_3, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_3, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_4, color=strain))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_4, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_5, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_5, color=cluster))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_5, color=sample))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_6, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_6, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_7, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_7, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_8, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_8, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_9, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_9, color=cluster))+geom_point()+theme_classic()

ggplot(U, aes(x=factor_1, y=factor_10, color=strain))+geom_point()+theme_classic()
ggplot(U, aes(x=factor_1, y=factor_10, color=cluster))+geom_point()+theme_classic()

head(V[order(V$factor_9, decreasing = T),], 20)
head(V[order(V$factor_5, decreasing = F),], 20)
head(V[order(V$factor_2, decreasing = F),], 20)

#### define another model and include sample information as a covariate as well

######## changing the cell_iD of U
sample_info = sapply(strsplit(U$cell_id, '_'), function(x) paste0(x[2], '_' ,x[3]))
sample_info[sample_info=='DA_M'] = 'DA_02'
sample_info[sample_info=='Lew_01'] = 'LEW_01'
sample_info[sample_info=='Lew_02'] = 'LEW_02'

table(sample_info)
pure_cell_id = sapply(strsplit(U$cell_id, '_'), function(x) paste0(x[length(x)]))
U$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
plot(U[,1],U[,2])


umap_df = data.frame(merged_samples[["umap"]]@cell.embeddings)
umap_df$cell_id = processed_df$cell_id
umap_df$cell_id_2 = processed_df$cell_id_2

sum(!U$cell_id_2 %in% umap_df$cell_id_2)
sum(!umap_df$cell_id_2 %in% U$cell_id_2)
head(umap_df$cell_id_2[!umap_df$cell_id_2 %in% U$cell_id_2])
merged_df = merge(U, umap_df, by.x='cell_id_2', by.y='cell_id_2')
dim(U)
dim(merged_df)
head(merged_df)

Var_to_check = 'factor_10'
ggplot(merged_df, aes(UMAP_1, UMAP_2, color=factor_10))+geom_point()+theme_classic()+
  scale_color_viridis(option = 'inferno', direction = -1)+ggtitle(Var_to_check)+
  theme(plot.title = element_text(size = 20))



