source('~/RatLiver/Codes/Functions.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(lme4) # load library
library(arm) # convenience functions for regression in R
library(dplyr)

preprocess_data <- function(a_dataset){
  
  a_dataset = data.frame(a_dataset)
  ### set the column names
  colnames(a_dataset) = a_dataset[1,]
  a_dataset = a_dataset[-1,]
  ### replacing the NA gene names with IDs
  tochange = which(is.na(a_dataset$AssociatedGeneName))
  a_dataset$AssociatedGeneName[tochange] = a_dataset$EnsemblGeneID[tochange]
  ## set the gene names
  rownames(a_dataset) = make.unique(a_dataset$AssociatedGeneName)
  ### create a seurat object
  a_dataset.s = CreateSeuratObject(counts = (a_dataset[,-c(1:4)]))
  ### add gene related meta data
  a_dataset.s[['RNA']]@meta.features$EnsemblGeneID = a_dataset$EnsemblGeneID
  a_dataset.s[['RNA']]@meta.features$EnsemblTranscriptID = a_dataset$EnsemblTranscriptID
  a_dataset.s[['RNA']]@meta.features$GeneLength = a_dataset$GeneLength
  
  return(a_dataset.s)
}


getResDF <- function(Model){
  df = data.frame(fitted=fitted(Model),
                  residuals=residuals(Model),
                  sample=concat_samples$sample,
                  PC1=Embeddings(concat_samples, 'pca')[,1])
  return(df)
}


data_files = list.files('~/scLMM/scLVM_Data/processed_Data/', pattern = '*singlecells_counts.txt', full.names = T)
data_list = lapply(data_files, read.table) 
names(data_list)  = c('G1', 'G2M', 'S')
lapply(data_list, dim)

seur_list = lapply(data_list, preprocess_data)


######################################################
############# Normalizing each sample individually ##########
# seur_norm_list <- lapply(X = seur_list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = seur_norm_list)
# length(features)
# 
# 
# concat_samples <- merge(seur_norm_list[[1]], c(seur_norm_list[[2]], seur_norm_list[[3]]), # 
#                         add.cell.ids = names(seur_list), 
#                         project = names(seur_list), 
#                         merge.data = TRUE)
# 
# VariableFeatures(concat_samples) = features
##################################################

############# Normalizing the merged samples ##########
concat_samples <- merge(seur_list[[1]], c(seur_list[[2]], seur_list[[3]]), # 
                        add.cell.ids = names(seur_list), 
                        project = names(seur_list), 
                        merge.data = TRUE)

#concat_samples <- NormalizeData(concat_samples) ## leads to dense cell populations
concat_samples <- SCTransform(concat_samples, variable.features.n = nrow(concat_samples))

concat_samples <- ScaleData(concat_samples)
#concat_samples = FindVariableFeatures(concat_samples, selection.method = "vst", nfeatures = 2000)

concat_samples <- RunPCA(concat_samples, verbose = FALSE, features=rownames(concat_samples))#features = features 
ElbowPlot(concat_samples)
concat_samples$sample = unlist(lapply(str_split(colnames(concat_samples), '_'), '[[', 1))

df_pca = data.frame(Embeddings(concat_samples, 'pca'), 
                    sample=concat_samples$sample)
ggplot(df_pca, aes(x=PC_1, y=PC_2, color=sample))+geom_point(alpha=0.8, size=2)+theme_classic()


PCA_loadings = data.frame(PC1_loading=Loadings(concat_samples, 'pca')[,1])
PCA_loadings$gene = rownames(PCA_loadings)
head(PCA_loadings)
dim(PCA_loadings)
##################################################

#saveRDS(df_pca, '~/scLMM/scLVM_Data/PCA_res.rds')
#saveRDS(concat_samples, '~/scLMM/scLVM_Data/cellCycleSamplesSeurObj.rds')



#################################################################### 
################# LMM model on cell-cycle ################# 
data = GetAssayData(concat_samples)

df = data.frame(sample=concat_samples$sample,
                PC1=Embeddings(concat_samples, 'pca')[,1],
                PC2=Embeddings(concat_samples, 'pca')[,2])

head(df)
hist(df$PC1)
ggplot(df, aes(x=PC1))+geom_density()+theme_classic()


############################################################
######### Binning the PC1 to make it categoriacal, so it could be used as a categorical variable ########## 
############################################################
# library(logiBin)
# b1 <- getBins(df = input.data, y='PC1', xVars=c('gene', 'sample'), minCr = 0.8, nCores = 2)
# 
# library(rbin)
# bins <- rbin_quantiles(input.data, response = gene, predictor = PC1, bins = 20)
# bins
# plot(bins)
# hist(input.data$gene)
# hist(input.data$PC1)
# 
# rbin_quantiles(mbank, response = y, predictor = day, bins = 10)
# hist(mbank$age)


library("optbin")
bin_thr = optbin(x=df$PC1, numbin=10, metric=c('se', 'mse'), is.sorted=FALSE, max.cache=2^31, na.rm=FALSE)
df$PC1_bin = assign.optbin(x=df$PC1, binspec=bin_thr, extend.upper=FALSE, by.value=FALSE)
plot(bin_thr)
hist(bin_thr)

df_pca$PC1_bin = as.character(df$PC1_bin)
ggplot(df_pca, aes(x=PC_1, y=PC_2, color=PC1_bin))+geom_point(alpha=0.8, size=2)+theme_classic()


Fit_lmer_model <- function(i){
  input.data = data.frame(gene=t(data)[,i], #'Itgal'
                          sample=as.character(as.factor(concat_samples$sample)),
                          PC1=Embeddings(concat_samples, 'pca')[,1],
                          numFeatures=concat_samples$nFeature_RNA,
                          librarySize=concat_samples$nCount_RNA, 
                          PC1_bin=df$PC1_bin)
  
  Model = lmer(gene ~ (1|PC1_bin), data=input.data, REML=T) 
  return(Model)
}
############################################################


Fit_lmer_model_2 <- function(i){
  input.data = data.frame(gene=t(data)[,i], #'Itgal'
                          sample=as.character(as.factor(concat_samples$sample)),
                          PC1=Embeddings(concat_samples, 'pca')[,1],
                          numFeatures=concat_samples$nFeature_RNA,
                          librarySize=concat_samples$nCount_RNA)
  
  Model = lmer(gene ~ (1|sample), data=input.data, REML=T) 
  return(Model)
}


start <- Sys.time()
models <- sapply(1:nrow(data), function(i) Fit_lmer_model_2(i), simplify = F)
end <- Sys.time()
total_time <- as.numeric(end - start, units = "mins") 
print(total_time)


#models = readRDS('~/scLMM/PCA_binned_lmer_models.rds')
#models = readRDS('~/scLMM/PCA_binned_lmer_models_allGenesPCA.rds')
models = readRDS('~/scLMM/sampleName_lmer_models_allGenesPCA.rds')


names(models) <- rownames(data)

is_model_converged <- unlist(lapply(models, function(x) isEmpty(unlist((x@optinfo$conv$lme4)))))
sum(!is_model_converged)*100/length(is_model_converged)

get_RandEff_Var <- function(model){
  sum = summary(model)
  RandomEffVar = data.frame(sum$varcor) # sdcor: SD - vcov: Variance
  return(RandomEffVar$vcov[1]) 
}

get_Residual_Var <- function(model){
  sum = summary(model)
  RandomEffVar = data.frame(sum$varcor) # sdcor: SD - vcov: Variance
  return(RandomEffVar$vcov[2]) 
}

RandEff_Var_list = lapply(models, get_RandEff_Var)
RandVar_df = data.frame(gene=names(RandEff_Var_list), Rand_Var=unlist(RandEff_Var_list))
RandVar_df_ord = RandVar_df[order(RandVar_df$Rand_Var, decreasing = T),]
rownames(RandVar_df_ord) = NULL
gridExtra::grid.table(head(RandVar_df_ord,20))
dev.off()

Residual_Var_list = lapply(models, get_Residual_Var)
Residual_Var_df = data.frame(gene=names(Residual_Var_list), Resid_Var=unlist(Residual_Var_list))
Residual_Var_df_ord = Residual_Var_df[order(Residual_Var_df$Resid_Var, decreasing = T),]
rownames(Residual_Var_df_ord) = NULL
gridExtra::grid.table(head(Residual_Var_df_ord,20))


Variance_df = cbind(RandVar_df, Residual_Var_df, is_model_converged)[,-1]
head(Variance_df)
ggplot(Variance_df, aes(Rand_Var, Resid_Var, color=Rand_Var))+geom_point()+theme_classic()+scale_color_viridis(option = 'magma',direction = 1)
ggplot(Variance_df, aes(Rand_Var, Resid_Var, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))

Variance_df = cbind(RandVar_df, Residual_Var_df)[,-1]
merged_df = merge(Variance_df, PCA_loadings, by.x='gene', by.y='gene')
dim(merged_df)
merged_df$PC1_loading_abs <- abs(merged_df$PC1_loading)
pheatmap(cor(merged_df[,-1]))

############# Rand_Var and PC loading
ggplot(merged_df, aes(Rand_Var, PC1_loading_abs, color=Resid_Var))+geom_point()+theme_classic()+scale_color_viridis(option = 'magma',direction = 1)
ggplot(merged_df, aes(Rand_Var, PC1_loading, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
ggplot(merged_df, aes(Rand_Var, PC1_loading_abs, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
############# Resid_Var and PC loading
ggplot(merged_df, aes(Resid_Var, PC1_loading, color=Rand_Var))+geom_point()+theme_classic()+scale_color_viridis(option = 'inferno',direction = 1)
ggplot(merged_df, aes(Resid_Var, PC1_loading, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
ggplot(merged_df, aes(Resid_Var, PC1_loading_abs, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))

#############
ggplot(merged_df, aes(Rand_Var, Resid_Var , color=PC1_loading_abs))+geom_point()+theme_classic()+scale_color_viridis(option = 'magma',direction = 1)


#### checking the model which had highest amount of explained variance by PC1
model_2test =  models[[RandVar_df_ord[1,]$gene]]
coef_test = coef(model_2test)
coef_test$sample

######################################## 
##### Running Coeff on the models ##### 
########################################
Models_coef_list = lapply(models, coef)
models_coef_df = Reduce(f = cbind, lapply(Models_coef_list, function(x) x$sample)) #PC1_bin
colnames(models_coef_df) = names(models)
rownames(models_coef_df) = paste0('PC_bin_',rownames(models_coef_df))
rownames(models_coef_df) = rownames(models_coef_df)

models_coef_df = t(models_coef_df)
gridExtra::grid.table(head(round(models_coef_df,4), 20))
dev.off()
pheatmap(cor(models_coef_df))

models_coef_df = data.frame(models_coef_df)
a_rand_coef_tab = models_coef_df[order(models_coef_df$S, decreasing = T),]
gridExtra::grid.table(head(round(a_rand_coef_tab,4), 20))

##############################################
df_3 = getResDF(Model)
ggplot(df_3, aes(x=fitted, y=residuals, color=sample))+
  geom_point(alpha=0.8, size=2.5)+theme_classic()


hist(df_3$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model))
coef(Model)
