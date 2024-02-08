source("~/scLMM/LMM-scRNAseq/lmmfit.R")
source("~/scLMM/LMM-scRNAseq/lmmfitSS.R")
source("~/scLMM/LMM-scRNAseq/lmmtest.R")
source("~/scLMM/LMM-scRNAseq/qqpvalue.R")

dirWork = '~/scLMM/LMM-scRNAseq/'
setwd(dirWork)
dir(dirWork)

library(Matrix)
library(MASS)
#BiocManager::install('maditr')
library(maditr)
library(data.table)
library(Seurat)

#############################################
######### Working with normalized Datasets 
#############################################

kidney <- readRDS('~/scLMM/LMM-scRNAseq/Data/Human_Kidney_data.rds')
dim(kidney)
dim(kidney)

kidney = kidney[rowSums(kidney)!=0,colSums(kidney)!=0]
Y_kidney = GetAssayData(kidney)
Y_kidney[100:120,100:150]
nGenes <- colSums(Y_kidney)
#Ylog_kidney = log2(t(Y_kidney) + 1)
#Ylog_kidney <- as.matrix(Ylog_kidney)
Ylog_kidney <- as.matrix(t(Y_kidney))
head(Ylog_kidney[1:10,1:10])


cell_type_list = names(table(kidney$Cell_Types_Broad))
split_gene_exp = sapply(1:length(cell_type_list), 
       function(i) GetAssayData(kidney)[,kidney$Cell_Types_Broad %in% cell_type_list[i]], simplify = F)

split_gene_exp_ave = sapply(1:length(split_gene_exp), 
                            function(i) data.frame(rowSums(split_gene_exp[[i]])/ncol(split_gene_exp[[i]])), simplify = F)
names(split_gene_exp_ave) = cell_type_list
lapply(split_gene_exp_ave, head)
gene_exp_ave_df = data.frame(Reduce(cbind, split_gene_exp_ave))
head(gene_exp_ave_df)
colnames(gene_exp_ave_df) = cell_type_list
pheatmap::pheatmap(cor(gene_exp_ave_df))
colnames(cor(gene_exp_ave_df))

kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c( "CCD-like" , "LOH-like"), "CCD-like", kidney$Cell_Types_Broad)
kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c("IC-A", "IC-B" ), "IC", kidney$Cell_Types_Broad)
kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c("DCT", "cTAL" ), "DCT", kidney$Cell_Types_Broad)
kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c('PC', 'CNT'), "PC", kidney$Cell_Types_Broad)
kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c('MNP', 'NK cell', 'B cell', 'T cell'), "lymphocyte", kidney$Cell_Types_Broad)
kidney$Cell = ifelse(kidney$Cell_Types_Broad %in% c('PEC', 'U1', 'U2', 'Endothelial', "Mesangial"), "endo", kidney$Cell_Types_Broad)



##design matrix for random effects
kidney$Cell[is.na(kidney$Cell)] = 'unknown'
kidney$sex[is.na(kidney$sex)] = 'unknown'

table(kidney$Cell)
sum(table(kidney$Cell_Types_Broad)) == ncol(kidney)
sum(table(kidney$sampleID)) == ncol(kidney)
Z_kidney = data.frame(ind=kidney$sampleID)
Z_kidney = model.matrix(~0+ind, Z_kidney) ### double check this with changjiang
#Z_pbmc = model.matrix(~0+ind, Z_pbmc) ### double check this with changjiang

head(Z_kidney,3)

##dimemsion of random effects
d_kidney <- ncol(Z_kidney)

##design matrix for fixed effects
X_kidney = data.frame(sex=kidney$sex,cell=kidney$Cell)
X_kidney = model.matrix(~ 0 + cell + sex:cell, X_kidney) #
X_kidney <- cbind(X_kidney, log_nGenes = log(nGenes))
unique(colnames(X_kidney))

#colnames(X_kidney)[1] = 'Intercept'
head(X_kidney)
class(X_kidney)

dim(X_kidney)
dim(Z_kidney)
dim(Ylog_kidney)


##############################
## **Fit LMM by lmmfit.**
#########################
SEED <- 57774
set.seed(SEED)
s0 <- c(runif(length(d_kidney)), 1)

t1 <- Sys.time()
fit <- NULL
fit <- lmmfit(Y = Ylog_kidney, X = X_kidney, Z = Z_kidney, 
              d = d_kidney, s0 = s0, max.iter = 50, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1) 

str(fit)
table(fit$niter)
sum(abs(fit$dlogL) > 1e-3)
dim(Ylog_kidney)


##############################
## Fit LMM by lmmfitSS.**
##############################

X = X_kidney; Y= Ylog_kidney; Z= Z_kidney; d=d_kidney
##input for lmmfitSS:
##- Compute XY, ZX, ZY, ZZ, XXinv, Ynorm, n
n <- nrow(X)
XY <- as.matrix(t(X)%*%Y)
ZX <- as.matrix(t(Z)%*%X)
ZY <- as.matrix(t(Z)%*%Y)
ZZ <- as.matrix(t(Z)%*%Z)
XXinv <- as.matrix(ginv(t(X)%*%X))
Ynorm <- colSums(Y*Y)

rm(X, Y, Z)

##Run lmmfitSS
set.seed(SEED)
s0 <- c(runif(length(d)), 1)

t1 <- Sys.time()
fitss <- lmmfitSS(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, 
                  Ynorm = Ynorm, n = n, d = d, 
                  s0 = s0, max.iter = 50, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1) 

#str(fitss)
##Compare lmmfit and lmmfitSS
range(fit$dlogL - fitss$dlogL) 
range(fit$niter - fitss$niter) 
range(fit$theta - fitss$theta) 
range(fit$coef - fitss$coef)
range(fit$cov - fitss$cov)   
##############################

##############################
##Test all coefficients.
##############################
##cluster k vs cluster 10 (intercept)
lmm <- lmmtest(fitss)
coef_df = data.frame(t(fitss$coef))
head(coef_df)

lmm_t = lmm[, grep("_t", colnames(lmm))]
lmm_pvalue = lmm[, grep("pvalue", colnames(lmm))]
head(lmm_t)
head(lmm_pvalue)

info_df = sapply(1:ncol(lmm_pvalue), function(i){
  df=data.frame(coef=coef_df[,i],
                tval=lmm_t[,i],
                pval=lmm_pvalue[,i])
  df$score = -log(df$pval+1e-16)*df$coef
  df=df[order(df$score, decreasing = T),]
  return(df)
}, simplify = F)
names(info_df) = colnames(lmm_pvalue)
lapply(info_df[10:18], head, 30)
lapply(info_df[16:30], head, 30)

names(info_df)
