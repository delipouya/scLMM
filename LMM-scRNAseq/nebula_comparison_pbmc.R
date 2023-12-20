library(Matrix)
library(MASS)
library(maditr)
library(data.table)
library(Seurat)
library(nebula)


pbmc <- readRDS("~/scLMM/LMM-scRNAseq//Data/PBMC_Lupus_Kang8vs8_data_counts.rds")
dim(pbmc)
pbmc = pbmc[rowSums(pbmc)!=0,colSums(pbmc)!=0]
dim(pbmc)
Y_pbmc = GetAssayData(pbmc)
nGenes <- colSums(Y_pbmc)
Ylog_pbmc = log2(t(Y_pbmc) + 1)
Ylog_pbmc <- as.matrix(Ylog_pbmc)

pbmc$cell[is.na(pbmc$cell)] = 'unknown'
table(pbmc$cell)

#################################
##design matrix for fixed effects
X_pbmc = data.frame(stim=pbmc$stim,cell=pbmc$cell)
X_pbmc = model.matrix(~ 1 + cell + stim, X_pbmc)
head(X_pbmc)

data_g = group_cell(count= as.matrix(Y_pbmc),id=pbmc$ind,pred=X_pbmc)
t1 <- Sys.time()
re = nebula(data_g$count,data_g$id,pred=data_g$pred)
t2 <- Sys.time()

print(paste0('nebula model I run times:', difftime(t2, t1)))
saveRDS(re, '~/scLMM/LMM-scRNAseq/NEBULA_modelI_stimPBMC.rds')


#################################
##design matrix for fixed effects

X_pbmc = data.frame(stim=pbmc$stim,cell=pbmc$cell)
X_pbmc = model.matrix(~ 1 + cell + stim:cell, X_pbmc)
head(X_pbmc)
num_cells_to_include = ncol(Y_pbmc)

data_g = group_cell(count= as.matrix(Y_pbmc),id=pbmc$ind,pred=X_pbmc)
t1 <- Sys.time()
re = nebula(data_g$count,data_g$id,pred=data_g$pred)
t2 <- Sys.time()

difftime(t2, t1)

print(paste0('nebula model II run times:', difftime(t2, t1)))
saveRDS(re, '~/scLMM/LMM-scRNAseq/NEBULA_modelII_stimPBMC.rds')

