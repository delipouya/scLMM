source('~/RatLiver/Codes/Functions.R')
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

### import the preprocessed input data 
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

#### generating information data frame based on the pre-processed/normalized data 
norm_df = data.frame(cell_ID=colnames(your_scRNAseq_data_object),
                     strain=as.factor(sapply(str_split(colnames(your_scRNAseq_data_object), '_'), '[[', 2)),
                     sample=as.factor(sapply(str_split(colnames(your_scRNAseq_data_object), '_'), 
                                              function(x){paste0(x[2],'_' ,x[3])})),
                     barcode=sapply(str_split(colnames(your_scRNAseq_data_object), '_'), function(x){x[length(x)]}),
                     clusters=sCVdata_list$res.0.2@Clusters)

norm_df$refined_sample = ifelse(norm_df$sample=='DA_M', 'DA_02',ifelse(norm_df$sample=='Lew_01', 'LEW_01', ifelse(norm_df$sample=='Lew_02', 'LEW_02', 'DA_01')))
table(norm_df$refined_sample, norm_df$sample)
norm_df$refined_ID = paste0(norm_df$refined_sample , '_', norm_df$barcode)
head(norm_df)

################################
set1_merged_samples_count = readRDS('~/scLMM/input_data_designMat/set1_merged_samples_count.rds')
table(sapply(str_split(colnames(set1_merged_samples_count), '_'), function(x){paste0(x[1],'_' ,x[2])}))


count_df = data.frame(cell_ID=colnames(set1_merged_samples_count),
                      strain=as.factor(sapply(str_split(colnames(set1_merged_samples_count), '_'), '[[', 1)),
                      barcode=sapply(str_split(colnames(set1_merged_samples_count), '_'), function(x){x[length(x)]}),
                      sample=as.character(sapply(str_split(colnames(set1_merged_samples_count), '_'), 
                                    function(x){paste0(x[1],'_' ,x[2])})))

table(count_df$sample)
count_df$sample[count_df$sample == 'DA_10WK'] = 'DA_02'
count_df$refined_ID = paste0(count_df$sample , '_', count_df$barcode)
head(count_df)


merged_info = merge(count_df, norm_df, by.x='refined_ID', by.y='refined_ID', all.x=F, sort=F) ### many cells were not matching well - need to be evaluated later
dim(merged_info)
head(merged_info)
merged_info2 = merged_info[,colnames(merged_info) %in% c('refined_ID', 'refined_sample', 'strain.x','clusters')]
sum(count_df$refined_ID != merged_info$refined_ID) ### the right order has been preserved
colnames(merged_info2) = c('ID','strain', 'cluster', 'sample')


#DA_01 DA_10WK  LEW_01  LEW_02 
#rat_DA_01   rat_DA_M rat_Lew_01 rat_Lew_02

design_mat <- model.matrix(~strain+sample+cluster, data = merged_info2)
#saveRDS(as.data.frame(design_mat), '~/scLMM/input_data_designMat/designMatrix_rat_set1_countData_refined.rds')
input.data = GetAssayData(set1_merged_samples_count)
#saveRDS(input.data, '~/scLMM/inputdata_rat_set1_countData.rds')

#####################  GLMM - Adding cluster as a random effect ##################### 
cluster_df = data.frame(clusters=sCVdata_list$res.0.2@Clusters, 
                        cell_ID=colnames(your_scRNAseq_data_object),
                        strain=sapply(str_split(colnames(your_scRNAseq_data_object), '_'), 
                                      function(x){paste0(x[1],'_' ,x[2],'_' ,x[3])}))

cluster_df$strain_2 = ifelse(cluster_df$strain=='rat_DA_01', 'DA_01', 
                  ifelse(cluster_df$strain=='rat_DA_M', 'DA_10WK', 
                         ifelse(cluster_df$strain=='rat_Lew_01', 'LEW_01', 'LEW_02')))

cluster_df$cell_ID_2 =  paste0(cluster_df$strain_2, '_', sapply(str_split(colnames(your_scRNAseq_data_object), '_'), '[[', 4))

cluster_df_sub = data.frame(ID=cluster_df$cell_ID_2, sample=cluster_df$strain_2, cluster=cluster_df$clusters)

head(count_df)
head(cluster_df)
test = merge(count_df, cluster_df_sub, by.x='cell_ID', by.y='ID', all.x=F) ### many cells were not matching well - need to be evaluated later
dim(test)

not_included = cluster_df_sub[!cluster_df_sub$ID %in% count_df$cell_ID,]
  
table(sapply(str_split(cluster_df_sub$ID, '_'), function(x){paste0(x[1],'_' , x[2])}))




##################### Data for LMM ##################### 
### generating the input data as a sparse matrix. Rows represent features (genes/phenotypes) and columns represent samples (cells)
### if the input needs to be converted to dense matrix, use as.matrix(input.data) to make the conversion
input.data = GetAssayData(your_scRNAseq_data_object)
saveRDS(input.data, '~/scLMM/inputdata_rat_set1.rds')

strain_info = sapply(str_split(colnames(your_scRNAseq_data_object), '_'), '[[', 2)
strain_info = sapply(str_split(colnames(set1_merged_samples_count), '_'), function(x){paste0(x[1],'_' ,x[2])})
table(strain_info)

Z.df = data.frame(cell_ID = colnames(your_scRNAseq_data_object),
                  strain = as.factor(strain_info)
                  )

Z.df$cluster = as.factor(paste0('cluster_', as.character(sCVdata_list$res.0.2@Clusters)))


table(cluster_df$strain)

###### new sample (set-2) strain markers werer flipped - fixing the bug #####
#your_scRNAseq_data_object$strain = 
### flipping the sample tag

table(Z.df$strain)
design_mat <- model.matrix(~strain+cluster, data = Z.df)
saveRDS(as.data.frame(design_mat), '~/scLMM/designMatrix_rat_set1.rds')





