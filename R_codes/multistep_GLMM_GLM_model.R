# Predict the gene expression using a mixed model (GLMM, LMM - could try both or even use LMM with the count data) y = Zu + e
# Extract the prediction of the model: d=Zu (shown below)
# Incorporate d in a GLM model -> y = Xb + d + e
# X include the strain information for each cell
# Evaluate b to see if it makes sense
# Does b give us the gene information here? I think so.
# Instead of the factor discovery step, we could
# Extract the residual matrix (e)
# Apply PCA to the residual and evaluate the latent variables
# Are they informative

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
library(reshape2)
get_head<- function(df) df[1:4, 1:4]

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

########################################################
######## importing the data and quality control ########
########################################################
##design matrix
design_mat = readRDS('~/scLMM/input_data_designMat/designMatrix_rat_set1_countData_refined.rds')
str(design_mat)
dim(design_mat) 
head(design_mat)

##raw counts
dat <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
sample_info = sapply(strsplit(colnames(dat), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
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

ncol(Y) == nrow(design_mat)
genelist = rownames(Y)

################################################################################
####### Calculating Zu using GLMM o be incorporated into a GLM model ########### 
################################################################################

fitted_val_df = data.frame(cell_id=colnames(Y))
family = "poisson"

# results_list = list()
#genelist_org = genelist
#genelist = genelist[1:10]

t1 <- Sys.time()
model_fit_errors = rep(0,length(genelist))

for (i in 1:length(genelist)){
  
  y = Y[genelist[i],]
  glmer_data = data.frame(y=y, rand=sample_info) #design_mat$strainLEW
  formula = 'y ~ 1 + (1 | rand)'
  
  glmm.fit <- glmer(formula = as.formula(formula), 
                     data = glmer_data, 
                     family = family)
  
  summary(glmm.fit)
  #results_list[[genelist[i]]] = glmm.fit
  fitted_val_df$gene_RE_fit = fitted(glmm.fit)
  colnames(fitted_val_df)[ncol(fitted_val_df)] = genelist[i]
  if(!is.null(glmm.fit@optinfo$conv$lme4$messages)) model_fit_errors[i]=1
}

t2 <- Sys.time()
print(model_fit_errors)
head(fitted_val_df[1:5])

#saveRDS(fitted_val_df, 'Zu_glmm_sampleRE_df.rds')
#saveRDS(model_fit_errors, 'error_geneindex_glmm_sampleRE.rds')

print(paste0('time passed: ', round(t2-t1,3), 's'))


################################################################################
####### Using the calculated d=Zu in a GLM with strain as covariate ########### 
################################################################################

fitted_val_df = readRDS('~/scLMM/Zu_glmm_sampleRE_df.rds')
model_fit_errors = readRDS('~/scLMM/error_geneindex_glmm_sampleRE.rds')

get_head(fitted_val_df)
dim(fitted_val_df)
sum(!is.numeric(colSums(fitted_val_df[,-1])))

paste0(round(sum(model_fit_errors)*100/length(model_fit_errors),2), '% of the models had a singular fit')
## 1461/12103 models had convergence error - 12.07% of the models had a singular fit

fitted_val_df_2 = fitted_val_df[,-1] ## removing the cell_id column
get_head(fitted_val_df_2)

res.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=ncol(Y)))
coef.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=3))
rownames(res.matrix) = rownames(Y)
rownames(coef.matrix) = rownames(Y)
colnames(coef.matrix) = c('Intercept','strainLEW','sample_eff')

strain_info = sapply(strsplit(colnames(dat), '_'), '[[', 1)

for (i in 1:length(genelist)){ #
  y = Y[genelist[i],]
  glm_data = data.frame(y=y, strain=strain_info, 
                        sample_eff=fitted_val_df_2[,i]) #design_mat$strainLEW
  glm.model = glm('y ~ strain + sample_eff', family="poisson", data=glm_data)
  #res = residuals(glm.model)
  #res.matrix[i,] = res
  coef.matrix[i,] = coef(glm.model)
}
#colnames(res.matrix) = colnames(Y)
#rownames(res.matrix) = rownames(Y)
dim(Y)

saveRDS(coef.matrix, '~/scLMM/coef_mat_glm_sample_strain_controled.rds')
#saveRDS(res.matrix, 'residual_mat_glm_sample_strain_controled.rds')


################################################################################
####### Using the calculated d=Zu in a GLM with strain as covariate ########### 
################################################################################

res.matrix = readRDS('~/scLMM/GLM_Results/multi_step/residual_mat_glm_sample_strain_controled.rds')
get_head(res.matrix)
summary(as.numeric(res.matrix[1000,])) ## checking the residual distribution in for a single gene/model

#pca.res = prcomp(res.matrix)

pca.res = readRDS('~/scLMM/GLM_Results/multi_step/residual_mat_pca_result.rds')
score_matrix = data.frame(pca.res$rotation)
loading_matrix = pca.res$x
get_head(score_matrix)
pca.res$sdev
screeplot(pca.res, npcs = 15)
residual_pca_df = data.frame(score_matrix[,1:20], cell_ID=row.names(score_matrix))

####### refining the cell ID information of the residual PCA results ####### 
sample_info = sapply(strsplit(residual_pca_df$cell_ID, '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
table(sample_info)
pure_cell_id = sapply(strsplit(residual_pca_df$cell_ID, '_'), function(x) paste0(x[length(x)]))
residual_pca_df$cell_id_2 = paste0(sample_info, '_', pure_cell_id)

##################################################################################################
###### importing the pre-analyzed sample information to use the cluster information ############## 
#################################################################################################
library(Seurat)
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
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

sum(!residual_pca_df$cell_id_2 %in% processed_df$cell_id_2)

head(residual_pca_df)
head(processed_df)

merged_df = merge(residual_pca_df, processed_df, by.x='cell_id_2', by.y='cell_id_2')
head(merged_df)
ggplot(merged_df, aes(PC1, PC2, color=cluster))+geom_point()+theme_bw() ### the first 8 PCs seem to be the most important
ggplot(merged_df, aes(PC1, PC2, color=strain))+geom_point()+theme_bw()
######################################################

plot(score_matrix[,1],score_matrix[,2])

### varimax PCA
num_PC = 20
res.scaled <- ScaleData(CreateSeuratObject(res.matrix))
res.scaled = GetAssayData(res.scaled)
get_head(res.scaled)
get_head(loading_matrix)
rot_data <- get_varimax_rotated(res.scaled, loading_matrix[,1:num_PC])

#rm(res.scaled)
rm(res.matrix)
gc()

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
pdf('~/scLMM/GLM_Results/multi_step/residual_varimax_onTop10PC_multiStep_GLM_GLMM.pdf', width = 14.5, height = 5)
for(i in 2:21){
  merged_df_var2$Var_to_check = merged_df_var[,i]
  p1=ggplot(merged_df_var2, aes(Varimax_1, Var_to_check, color=cluster))+
    geom_point()+theme_bw()+ylab(colnames(merged_df_var)[i]) ### the first 8 PCs seem to be the most important
  p2=ggplot(merged_df_var2, aes(Varimax_1, Var_to_check, color=strain))+
    geom_point()+theme_bw()+ylab(colnames(merged_df_var)[i])
  
  p3=ggplot(merged_df_var2, aes(UMAP_1, UMAP_2, color=Var_to_check))+
    geom_point()+theme_bw()+
    scale_color_viridis(option='inferno', direction=-1)
  gridExtra::grid.arrange(p1, p2,p3 ,ncol = 3)
}
dev.off()




####### evaluating the effect of nnumber of PCs chosen to run varimax on the residual matrix
num_PC = 5
rot_data <- get_varimax_rotated(res.scaled_matrix, loading_matrix[,1:num_PC])
rotated_scores_5 <- data.frame(rot_data$rotScores)
colnames(rotated_scores_5) = paste0('Varimax_', 1:ncol(rotated_scores_5), '.pc',num_PC )
head(rotated_scores_5)

num_PC = 10
rot_data <- get_varimax_rotated(res.scaled_matrix, loading_matrix[,1:num_PC])
rotated_scores_10 <- data.frame(rot_data$rotScores)
colnames(rotated_scores_10) = paste0('Varimax_', 1:ncol(rotated_scores_10), '.pc',num_PC )
head(rotated_scores_10)

num_PC = 20
rot_data <- get_varimax_rotated(res.scaled_matrix, loading_matrix[,1:num_PC])
rotated_scores_20 <- data.frame(rot_data$rotScores)
colnames(rotated_scores_20) =  paste0('Varimax_', 1:ncol(rotated_scores_20), '.pc',num_PC )
head(rotated_scores_20)


rotated_scores_merged = cbind(rotated_scores_5, rotated_scores_10)
pheatmap::pheatmap(cor(rotated_scores_merged)[1:5,6:ncol(rotated_scores_merged)])

rotated_scores_merged = cbind(rotated_scores_5, rotated_scores_20)
pheatmap::pheatmap(cor(rotated_scores_merged)[1:5,6:ncol(rotated_scores_merged)])

rotated_scores_merged = cbind(rotated_scores_10, rotated_scores_20)
pheatmap::pheatmap(cor(rotated_scores_merged)[1:10,11:ncol(rotated_scores_merged)])


########### comparing the varimax coordinates of the var on residuals with the var results on gene expression matrix
varimax_res_set1 <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1), '.orig')
score_set1$cell_id_2 = processed_df$cell_id_2


rotated_scores_5_tmp = rotated_scores_5
rotated_scores_5_tmp$cell_id = rownames(rotated_scores_5_tmp)
rotated_scores_merged = merge(rotated_scores_5_tmp, score_set1, by.x='cell_id', by.y='cell_id_2')
pheatmap::pheatmap(cor(rotated_scores_merged[,-1])[1:5,6:(ncol(rotated_scores_merged)-1)])


rotated_scores_10_tmp = rotated_scores_10
rotated_scores_10_tmp$cell_id = rownames(rotated_scores_10_tmp)
rotated_scores_merged = merge(rotated_scores_10_tmp, score_set1, by.x='cell_id', by.y='cell_id_2')
pheatmap::pheatmap(cor(rotated_scores_merged[,-1])[1:10,11:(ncol(rotated_scores_merged)-1)])

rotated_scores_20_tmp = rotated_scores_20
rotated_scores_20_tmp$cell_id = rownames(rotated_scores_20_tmp)
rotated_scores_merged = merge(rotated_scores_20_tmp, score_set1, by.x='cell_id', by.y='cell_id_2')
pheatmap::pheatmap(cor(rotated_scores_merged[,-1])[1:20,21:(ncol(rotated_scores_merged)-1)])


###  Comparing the 

coef.matrix <- readRDS('~/scLMM/coef_mat_glm_sample_strain_controled.rds')
head(coef.matrix)
head(coef.matrix[order(coef.matrix$strainLEW, decreasing = T),],20) ### LEW specific genes
head(coef.matrix[order(coef.matrix$strainLEW, decreasing = F),],20) ### DA specific genes
coef.matrix$genes = rownames(coef.matrix) 

var15 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC15_loadings.rnk')
colnames(var15) = c('gene', 'var15_load')
var5 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC5_loadings.rnk')
colnames(var5) = c('gene', 'var5_load')

merged_var15 = merge(coef.matrix, var15, by.x = 'genes', by.y='gene')
merged_varimax = merge(merged_var15, var5, by.x = 'genes', by.y='gene')
head(merged_varimax)
merged_varimax$convergence = as.character(merged_varimax$convergence)

ggplot(merged_varimax, aes(strainLEW, var15_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with GLM LEW strain coefficient estimates')+
  geom_text(aes(label=ifelse(var15_load>0.1 | var15_load<(-0.1)|strainLEW>10|strainLEW<(-10) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)

ggplot(merged_varimax, aes(strainLEW, var5_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with GLM LEW strain coefficient estimates')+
  geom_text(aes(label=ifelse(var5_load>0.1 | var5_load<(-0.1)|strainLEW>10|strainLEW<(-10) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)

ggplot(merged_varimax, aes(sample_eff, var5_load))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('Varimax 5 correlation with GLM sample-effect coefficient estimates\nRemoved 1278 rows containing missing values ')+
  geom_text(aes(label=ifelse(var5_load>0.1 | var5_load<(-0.1) ,as.character(genes),'')),
            hjust=0,vjust=0,size=3)



head(var15)
dim(var15)
var15[var15$gene == 'Faslg',]
var5[var5$gene == 'Faslg',]



