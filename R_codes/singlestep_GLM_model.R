####################################
###### in this script, poisson/NB GLM and LM models are trained on the rat liver dataset

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
library(lme4)
library(stringr)
library(MASS)
library(reshape2)

get_head<- function(df) df[1:4, 1:4]

########################################################
######## importing the data and quality control ########
########################################################
##raw counts
dat <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
sample_info = sapply(strsplit(colnames(dat), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
strain_info = sapply(strsplit(colnames(dat), '_'), '[[', 1)
table(sample_info)

dim(dat) #[1] 32883 23036
get_head(dat)
ncells <- rowSums(dat > 0)
hist(ncells)
hist(log2(ncells))
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))

genes_2keep = ncells >= 20
Y <- as.matrix(dat[genes_2keep, ])
genelist = rownames(Y)
rm(dat)
gc()


########################################################
############ constructing the design matrix ###########
########################################################
library(mltools)
library(data.table)
library('dbarts')


strain_design <- data.frame(one_hot(as.data.table(as.factor(strain_info))))
sample_design <- data.frame(one_hot(as.data.table(as.factor(sample_info))))
design = cbind(strain_design, sample_design)
colnames(design) = gsub('V1_', '', colnames(design))
design = design[,colnames(design) %in% c("LEW", "DA_01", "LEW_01")]
head(design)
design =  makeModelMatrixFromDataFrame(design)


###################################################################################
############ Fitting the model and saving the residuals and coefficients ##########
###################################################################################
res.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=ncol(Y)))
coef.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=4)) # 1+1(strain)+2(sample)
rownames(res.matrix) = rownames(Y)
rownames(coef.matrix) = rownames(Y)
colnames(coef.matrix) = c('Intercept','LEW_strain', 'DA_01_sample', 'LEW_01_sample')
family = 'poisson'

get_head(Y)

for (i in 1:length(genelist)){ #
  
  y = log(Y[genelist[i],]+1e-100) #LM data - log normalization
  #y = Y[genelist[i],] ## GLM data
  
  #model = glm(y ~ design + 1, family=family) # GLM poisson link function
  #model = glm.nb(y ~ design + 1) # GLM NB link function
  model = lm(y ~ design + 1) # LM model
  
  res = residuals(model)
  res.matrix[i,] = res
  coef.matrix[i,] = coef(model)
  
}

colnames(res.matrix) = colnames(Y)
rownames(res.matrix) = rownames(Y)
get_head(res.matrix)
head(coef.matrix)



######### The single-step process #########

#saveRDS(coef.matrix, '~/scLMM/GLM_Results/single_step/coef_mat_glm_sampleDirect_strain_controled_poisson.rds') # GLM (strain+sample) - poisson 
#saveRDS(res.matrix, '~/scLMM/GLM_Results/single_step/residual_mat_glm_sampleDirect_strain_controled_poisson.rds') # GLM (strain+sample) - poisson

#saveRDS(coef.matrix, '~/scLMM/GLM_Results/single_step/coef_mat_glm_sampleDirect_strain_controled_NB.rds') # GLM (strain+sample) - NB  
#saveRDS(res.matrix, '~/scLMM/GLM_Results/single_step/residual_mat_glm_sampleDirect_strain_controled_NB.rds') # GLM (strain+sample) - NB 

#saveRDS(coef.matrix, '~/scLMM/GLM_Results/single_step/coef_mat_lm_sample_strain_logNorm.rds') # LM (strain+sample) - log norm data  
#saveRDS(res.matrix, '~/scLMM/GLM_Results/single_step/residual_mat_lm_sample_strain_logNorm.rds') # LM (strain+sample) - log norm data 


################################################################################
####### Applying factor analysis to the resulting residuals  ########### 
################################################################################
res.matrix = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_glm_sampleDirect_strain_controled_poisson.rds') # GLM (strain+strain) - poisson
#res.matrix = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_glm_sampleDirect_strain_controled_NB.rds') # GLM (strain+strain) - NB
res.matrix = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_lm_sample_strain_logNorm.rds') # LM (strain+sample) - log norm data 

get_head<- function(df) df[1:4, 1:4]
get_head(res.matrix)
summary(as.numeric(res.matrix[1000,])) ## checking the residual distribution in for a single gene/model



pca.res = prcomp(res.matrix)
#saveRDS(pca.res, '~/scLMM/GLM_Results/single_step/residual_mat_pca_glm_sampleDirect_strain_controled_poisson.rds')
#saveRDS(pca.res, '~/scLMM/GLM_Results/single_step/residual_mat_pca_glm_sampleDirect_strain_controled_NB.rds')
#saveRDS(pca.res, '~/scLMM/GLM_Results/single_step/residual_mat_pca_lm_sample_strain_logNorm.rds')

###### single-step poisson


#pca.res = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_pca_glm_sampleDirect_strain_controled_poisson.rds')
#pca.res = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_pca_glm_sampleDirect_strain_controled_NB.rds')
pca.res = readRDS('~/scLMM/GLM_Results/single_step/residual_mat_pca_lm_sample_strain_logNorm.rds')


score_matrix = data.frame(pca.res$rotation)
loading_matrix = pca.res$x
get_head(score_matrix)

num_PCs = 20 #### re-run the varimax with this value
### varimax PCA

res.scaled <- ScaleData(CreateSeuratObject(res.matrix))
res.scaled = GetAssayData(res.scaled)
get_head(res.scaled)
get_head(loading_matrix)
rot_data <- get_varimax_rotated(res.scaled, loading_matrix[,1:num_PCs])

rm(res.scaled)
rm(res.matrix)
gc()
##################################################################################################
###### importing the pre-analyzed sample information to use the cluster information ############## 
#################################################################################################

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
processed_df$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
head(processed_df)

###########################################################################################
################# evaluating the non-rotated PCA results ##################################
####### refining the cell ID information of the residual PCA results ####### 
screeplot(pca.res, npcs = 20)
residual_pca_df = data.frame(score_matrix[,1:num_PCs], cell_ID=row.names(score_matrix))

sample_info = sapply(strsplit(residual_pca_df$cell_ID, '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
table(sample_info)
pure_cell_id = sapply(strsplit(residual_pca_df$cell_ID, '_'), function(x) paste0(x[length(x)]))
residual_pca_df$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
plot(score_matrix[,1],score_matrix[,2])

sum(!residual_pca_df$cell_id_2 %in% processed_df$cell_id_2)

head(residual_pca_df)
head(processed_df)

merged_df = merge(residual_pca_df, processed_df, by.x='cell_id_2', by.y='cell_id_2')
head(merged_df)
ggplot(merged_df, aes(PC1, PC4, color=cluster))+geom_point()+theme_bw() ### the first 8 PCs seem to be the most important
ggplot(merged_df, aes(PC1, PC2, color=strain))+geom_point()+theme_bw()
###########################################################################################


#######################################################################################
################# evaluating the rotated PCA results ##################################

rotated_loadings <- rot_data$rotLoadings
rotated_scores <- data.frame(rot_data$rotScores)
colnames(rotated_scores) = paste0('Varimax_', 1:ncol(rotated_scores))

sum(rownames(rotated_scores) != residual_pca_df$cell_ID)
rotated_scores$cell_id = residual_pca_df$cell_ID
rotated_scores$cell_id_2 = residual_pca_df$cell_id_2

merged_df_var = merge(rotated_scores, processed_df, by.x='cell_id_2', by.y='cell_id_2')
head(merged_df_var)

umap_df = data.frame(merged_samples[["umap"]]@cell.embeddings)
umap_df$cell_id = processed_df$cell_id
umap_df$cell_id_2 = processed_df$cell_id_2


merged_df_var2 = merge(merged_df_var, umap_df, by.x='cell_id_2', by.y='cell_id_2')
head(merged_df_var2)

library(viridis)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_glm_sampleDirect_strain_controled_poisson.pdf', width = 14.5, height = 5)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_glm_sampleDirect_strain_controled_NB.pdf', width = 14.5, height = 5)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_lm_sample_strain_logNorm.pdf', width = 14.5, height = 5)

for(i in 2:21){
  merged_df_var2$Var_to_check = merged_df_var[,i]
  p1=ggplot(merged_df_var2, aes(Varimax_1, Var_to_check, color=cluster))+
    geom_point()+theme_bw()+ylab(colnames(merged_df_var)[i]) ### the first 8 PCs seem to be the most important
  p2=ggplot(merged_df_var2, aes(Varimax_1, Var_to_check, color=strain))+
    geom_point()+theme_bw()+ylab(colnames(merged_df_var)[i])
  
  p3=ggplot(merged_df_var2, aes(UMAP_1, UMAP_2, color=Var_to_check))+
    geom_point()+theme_bw()+
    scale_color_viridis(option='inferno', direction=-1)
  gridExtra::grid.arrange(p1, p2, p3 ,ncol = 3)
}
dev.off()


########### comparing the varimax coordinates of the var on residuals with the var results on gene expression matrix
varimax_res_set1 <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1), '.orig')
score_set1$cell_id_2 = processed_df$cell_id_2

rotated_scores_20_tmp = rotated_scores
rotated_scores_20_tmp$cell_id = rownames(rotated_scores_20_tmp)
rotated_scores_merged = merge(rotated_scores_20_tmp, score_set1, by.x='cell_id', by.y='cell_id_2')
rotated_scores_merged = rotated_scores_merged[,!colnames(rotated_scores_merged) %in% c('cell_id','cell_id_2')]
head(rotated_scores_merged)
pheatmap::pheatmap(cor(rotated_scores_merged)[1:20,21:ncol(rotated_scores_merged)])


###  Comparing the 
coef.matrix = readRDS('~/scLMM/GLM_Results/single_step/coef_mat_glm_sampleDirect_strain_controled_poisson.rds') # GLM (strain+strain) - poisson 
coef.matrix = readRDS('~/scLMM/GLM_Results/single_step/coef_mat_glm_sampleDirect_strain_controled_NB.rds') # GLM (strain+strain) - NB 
coef.matrix <- readRDS('~/scLMM/GLM_Results/single_step/coef_mat_lm_sample_strain_logNorm.rds')

head(coef.matrix)
head(coef.matrix[order(coef.matrix$LEW_strain, decreasing = T),],20) ### LEW specific genes
head(coef.matrix[order(coef.matrix$LEW_strain, decreasing = F),],20) ### DA specific genes
coef.matrix$genes = rownames(coef.matrix) 

var15 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC15_loadings.rnk')
colnames(var15) = c('gene', 'var15_load')
var5 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC5_loadings.rnk')
colnames(var5) = c('gene', 'var5_load')

merged_var15 = merge(coef.matrix, var15, by.x = 'genes', by.y='gene')
merged_varimax = merge(merged_var15, var5, by.x = 'genes', by.y='gene')
head(merged_varimax)

ggplot(merged_varimax, aes(LEW_strain, var15_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with GLM LEW strain coefficient estimates')+
  geom_text(aes(label=ifelse(var15_load>0.1 | var15_load<(-0.1)|LEW_strain>45|LEW_strain<(-45) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)

ggplot(merged_varimax, aes(LEW_strain, var5_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with LM LEW strain coefficient estimates')+
  geom_text(aes(label=ifelse(var5_load>0.1 | var5_load<(-0.1)|LEW_strain>45|LEW_strain<(-45) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)

ggplot(merged_varimax, aes(DA_01_sample, var5_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('Varimax 5 correlation with LM sample-effect coefficient estimates\nRemoved 1278 rows containing missing values ')+
  geom_text(aes(label=ifelse(var5_load>0.1 | var5_load<(-0.1) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)



head(var15)
dim(var15)
var15[var15$gene == 'Faslg',]
var5[var5$gene == 'Faslg',]
