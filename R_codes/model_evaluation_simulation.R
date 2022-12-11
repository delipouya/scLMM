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

######### ######### ######### 
###### Encoding batch information within the design matrix (controling) 
###### Expect to identify group information in the residuals

######### defining the design matrix #########

batch_design <- data.frame(one_hot(as.data.table(as.factor(meta_data$Batch))))
group_design <- data.frame(one_hot(as.data.table(as.factor(meta_data$Group))))
design = cbind(batch_design, group_design)
colnames(design) = gsub('V1_', '', colnames(design))

indices_toInclude = c(1)
design_cols = c("Batch1", "Batch2", "Group1", "Group2", "Group3")
design = data.frame(design[,colnames(design) %in% design_cols[indices_toInclude]])
colnames(design) = design_cols[indices_toInclude]
head(design)
design =  makeModelMatrixFromDataFrame(data.frame(design))


###################################################################################
############ Fitting the model and saving the residuals and coefficients ##########
###################################################################################
Y = counts(sim)
res.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=ncol(Y)))
coef.matrix = data.frame(matrix(NA,nrow=nrow(Y), ncol=(ncol(design)+1))) # 1+1(strain)+2(sample)
rownames(res.matrix) = rownames(Y)
rownames(coef.matrix) = rownames(Y)
colnames(coef.matrix) = c('Intercept',colnames(design))
family = 'poisson'

get_head(Y)
genelist = rownames(Y)[1:50]

for (i in 1:length(genelist)){ 
  
  y = Y[genelist[i],] ## GLM data
  
  if(MODEL == 'POIS_GLM' ) {
    model = glm(y ~ design + 1, family=family) # GLM poisson link function
    
  } else if(MODEL == 'NB_GLM' ) {
    model = glm.nb(y ~ design + 1) # GLM NB link function
    
    } else if(MODEL == 'LOG_LM'){
    y = log(Y[genelist[i],]+1e-100) #LM data - log normalization
    model = lm(y ~ design + 1) # LM model
    
    }
  
  res = residuals(model)
  res.matrix[i,] = res
  coef.matrix[i,] = coef(model)
  
}

colnames(res.matrix) = colnames(Y)
rownames(res.matrix) = rownames(Y)s
get_head(res.matrix)

#### to evaluate the performance of GLMM-GLM models, we need to have at least 3 sources of variation 
# a technical one to account for using mixed model (the fitted values will be used to control for in the GLM)
# two biological ones - one will be used as a covariate in the GLM - like strain
# the second biological covariate will be designed to be captured within the residuals

## can we look into the residuals of mixed models? how reliable are they?

if (MODEL == 'POIS_GLMM'){
  model <- glmer(formula = y ~ 1 + (1 | design), family = family)
  summary(glmm.fit)
  #results_list[[genelist[i]]] = glmm.fit
  fitted_val_df$gene_RE_fit = fitted(glmm.fit)
  colnames(fitted_val_df)[ncol(fitted_val_df)] = genelist[i]
  if(!is.null(glmm.fit@optinfo$conv$lme4$messages)) model_fit_errors[i]=1
}

################## Saving and importing data ############
#saveRDS(coef.matrix, paste0('~/scLMM/GLM_Results/single_step/coef_mat_simulation_',MODEL,'.rds')) 
#saveRDS(res.matrix, paste0('~/scLMM/GLM_Results/single_step/residual_mat_simulation_',MODEL,'.rds')) 

coef.matrix = readRDS(paste0('~/scLMM/GLM_Results/single_step/coef_mat_simulation_',MODEL,'.rds'))
res.matrix = readRDS(paste0('~/scLMM/GLM_Results/single_step/residual_mat_simulation_',MODEL,'.rds')) 
  
get_head(res.matrix)
head(coef.matrix,20)

#pca.res = prcomp(res.matrix)
# saveRDS(pca.res, paste0('~/scLMM/GLM_Results/single_step/residual_PCA_simulation_',MODEL,'.rds')) 
pca.res = readRDS(paste0('~/scLMM/GLM_Results/single_step/residual_PCA_simulation_',MODEL,'.rds')) 

################## PCA score and loadings ##################
score_matrix = data.frame(pca.res$rotation)
loading_matrix = pca.res$x
get_head(score_matrix)
num_PCs = 10 #### re-run the varimax with this value
sum(rownames(score_matrix) != meta_data$Cell)
score_matrix_2 = cbind(score_matrix[,1:num_PCs], meta_data)
sum(rownames(loading_matrix) != gene_info$Gene)

#######################################################################################
######## evaluating the reason for getting many NA values for the NB GLM model ########
### TODO: unresoved ??
is_na_row = rowSums(is.na(res.matrix)) > 0
sum(is_na_row)
res.matrix_na_row = res.matrix[is_na_row,]
dim(res.matrix_na_row)
dim(res.matrix)

sumExp_NA_genes = rowSums(Y[is_na_row,])
sumExp_non_NA_genes = rowSums(Y[! is_na_row,])
summary(sumExp_NA_genes)
summary(sumExp_non_NA_genes)

head(coef.matrix[rowSums(is.na(coef.matrix)) > 0,])
########################################################################

### varimax PCA

res.scaled <- ScaleData(CreateSeuratObject(res.matrix))
res.scaled = GetAssayData(res.scaled)
get_head(res.scaled)
get_head(loading_matrix)
rot_data <- get_varimax_rotated(res.scaled, loading_matrix[,1:num_PCs])

###########################################################################################
################# visualizing the non-rotated PCA results ##################################
####### refining the cell ID information of the residual PCA results ####### 
screeplot(pca.res, npcs = 50)
residual_pca_df = data.frame(score_matrix[,1:num_PCs], cell_ID=row.names(score_matrix))

sum(residual_pca_df$cell_ID != meta_data$Cell)
residual_pca_df = cbind(residual_pca_df, meta_data)
head(residual_pca_df)

ggplot(residual_pca_df, aes(PC1, PC2, color=Batch))+geom_point()+theme_bw() 
ggplot(residual_pca_df, aes(PC1, PC2, color=Group))+geom_point()+theme_bw() 
ggplot(residual_pca_df, aes(PC1, PC3, color=Group))+geom_point()+theme_bw() 
ggplot(residual_pca_df, aes(PC1, PC6, color=Group))+geom_point()+theme_bw() 
ggplot(residual_pca_df, aes(PC1, PC8, color=Group))+geom_point()+theme_bw() 


ggplot(residual_pca_df, aes(PC1, PC2, color=ExpLibSize))+geom_point()+theme_bw()
ggplot(residual_pca_df, aes(PC1, PC2, color=sizeFactor))+geom_point()+theme_bw()
###########################################################################################



############### varimax PCA score annd loadings ##################
rotated_loadings <- rot_data$rotLoadings
rotated_scores <- data.frame(rot_data$rotScores)
colnames(rotated_scores) = paste0('Varimax_', 1:ncol(rotated_scores))

sum(rownames(rotated_scores) != residual_pca_df$cell_ID)
rotated_scores = cbind(rotated_scores, meta_data)
head(rotated_scores)

umap_df = data.frame(reducedDim(sim,type='UMAP'))
sum(rownames(umap_df) != rownames(rotated_scores))
colnames(umap_df) = c('UMAP_1', 'UMAP_2')  
rotated_scores = cbind(rotated_scores, umap_df)
head(rotated_scores)
######################################################

### evaluating the matching based on PCA factors
loading_matrix_2 = cbind(loading_matrix[,1:num_PCs], gene_info)

### evaluating the matching based on Varimax-PCA factors
loading_matrix_2 = cbind(rotated_loadings[,1:num_PCs], gene_info)


loading_matrix_2 = loading_matrix_2[,!colnames(loading_matrix_2) %in% c('Gene','BaseGeneMean','OutlierFactor')]
head(loading_matrix_2)
cor_matrix_loading = cor(loading_matrix_2)
PC_colnames = paste0('PC', 1:num_PCs)
cor_matrix_loading_sub = cor_matrix_loading[!colnames(cor_matrix_loading) %in% PC_colnames,colnames(cor_matrix_loading) %in% PC_colnames]

pheatmap(cor_matrix_loading_sub,color=colorRampPalette(c("navy", "white", "red"))(50), display_numbers = TRUE)

PCs_tocheck = colSums(abs(cor_matrix_loading_sub))>median(colSums(abs(cor_matrix_loading_sub)))
columns_tocheck = c(PC_colnames[PCs_tocheck],colnames(cor_matrix_loading)[!colnames(cor_matrix_loading) %in% PC_colnames])

panel.points<-function(x,y) points(x,y,cex=0.2)
pairs(loading_matrix_2[,colnames(loading_matrix_2) %in% columns_tocheck], 
      lower.panel=panel.points, upper.panel=panel.points)


get_max_cor <- function(covariate_name) { max(abs(cor_matrix_loading_sub[covariate_name,]))}
score = (get_max_cor('DEFacGroup1')* get_max_cor('DEFacGroup2') * get_max_cor('DEFacGroup3'))/(get_max_cor('BatchFacBatch1')*get_max_cor('BatchFacBatch2'))
print(score)

#######################################################################################
################# Visualizing the rotated PCA results ##################################

library(viridis)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_glm_sampleDirect_strain_controled_poisson.pdf', width = 14.5, height = 5)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_glm_sampleDirect_strain_controled_NB.pdf', width = 14.5, height = 5)
#pdf('~/scLMM/GLM_Results/single_step/residual_varimaxTop20PC_lm_sample_strain_logNorm.pdf', width = 14.5, height = 5)

ggplot(rotated_scores, aes(Varimax_1, Varimax_7, color=Group))+geom_point()+theme_bw() 
ggplot(rotated_scores, aes(Varimax_1, Varimax_10, color=Group))+geom_point()+theme_bw() 

for(i in 2:21){
  rotated_scores$Var_to_check = rotated_scores[,i]
  p1=ggplot(rotated_scores, aes(Varimax_1, Var_to_check, color=Batch))+
    geom_point()+theme_bw()+ylab(colnames(rotated_scores)[i]) ### the first 8 PCs seem to be the most important
  p2=ggplot(rotated_scores, aes(Varimax_1, Var_to_check, color=Group))+
    geom_point()+theme_bw()+ylab(colnames(rotated_scores)[i])

  
  p3=ggplot(rotated_scores, aes(UMAP_1, UMAP_2, color=Var_to_check))+
    geom_point()+theme_bw()+
    scale_color_viridis(option='inferno', direction=-1)+ggtitle(colnames(rotated_scores)[i])
  gridExtra::grid.arrange(p1, p2, p3 ,ncol = 3)
}
dev.off()


###### checking the batch coefficient estimation with 
head(coef.matrix[order(coef.matrix$Batch1, decreasing = T),],20) ### LEW specific genes
coef.matrix$genes = rownames(coef.matrix) 
sum(coef.matrix$genes != gene_info$Gene)
coef.matrix = cbind(coef.matrix, gene_info)
head(coef.matrix)
ggplot(coef.matrix, aes(Batch1, BatchFacBatch1))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle(paste0('Batch1 Est-coefficient correlation with true value - ', MODEL))
  #geom_text(aes(label=ifelse(var15_load>0.1 | var15_load<(-0.1)|LEW_strain>45|LEW_strain<(-45) ,as.character(genes),'')),hjust=0,vjust=0,size=3)
ggplot(coef.matrix, aes(Batch1, BatchFacBatch2))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle(paste0('Batch1 Est-coefficient correlation with true value - ', MODEL))

cor(coef.matrix$Batch1, coef.matrix$BatchFacBatch1)
########

rotated_loadings = as.data.frame(rotated_loadings)
sum(rownames(rotated_loadings) != gene_info$Gene)
loading_matrix2 = cbind(rotated_loadings[,1:10], gene_info)
head(loading_matrix2)
#### varimax PC1 loading
ggplot(loading_matrix2, aes(PC1, DEFacGroup1))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC1 loading')+
  theme_classic()+ggtitle(paste0('Group-1 Est-coefficient correlation with true value - ', MODEL))
ggplot(loading_matrix2, aes(PC1, DEFacGroup2))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC1 loading')+
  theme_classic()+ggtitle(paste0('Group-2 Est-coefficient correlation with true value - ', MODEL))
ggplot(loading_matrix2, aes(PC1, DEFacGroup3))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC1 loading')+
  theme_classic()+ggtitle(paste0('Group-3 Est-coefficient correlation with true value - ', MODEL))

#### varimax PC2 loading
ggplot(loading_matrix2, aes(PC2, DEFacGroup1))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC2 loading')+
  theme_classic()+ggtitle(paste0('Group-1 Est-coefficient correlation with true value - ', MODEL))
ggplot(loading_matrix2, aes(PC2, DEFacGroup2))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC2 loading')+
  theme_classic()+ggtitle(paste0('Group-2 Est-coefficient correlation with true value - ', MODEL))
ggplot(loading_matrix2, aes(PC2, DEFacGroup3))+geom_point(size=1.5, alpha=0.5)+xlab('Varimax PC2 loading')+
  theme_classic()+ggtitle(paste0('Group-3 Est-coefficient correlation with true value - ', MODEL))



