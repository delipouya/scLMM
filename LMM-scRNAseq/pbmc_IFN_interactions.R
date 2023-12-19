source("~/scLMM/LMM-scRNAseq/lmmfit.R")
source("~/scLMM/LMM-scRNAseq/lmmfitSS.R")
source("~/scLMM/LMM-scRNAseq/lmmtest.R")
source("~/scLMM/LMM-scRNAseq/qqpvalue.R")
library(Matrix)
library(MASS)
#BiocManager::install('maditr')
library(maditr)
library(data.table)
library(Seurat)

##########################################
########## how to define the design matrix
##########################################
n <- 12
set.seed(2012)
cell <- as.factor(sample(1:3, size = n, replace = T))
trt <- as.factor(sample(1:2, size = n, replace = T))
design <- model.matrix(~ 0 + cell + trt:cell)
design[1:3,]
##########################################


dirWork = '~/scLMM/LMM-scRNAseq/'
setwd(dirWork)
dir(dirWork)



#############################################
######### Working with normalized Datasets 
#############################################

##### PBMC dataset #####
#pbmc <- readRDS('LMM-scRNAseq/Data/PBMC_Lupus_Kang8vs8_data_norm.rds')
pbmc <- readRDS("~/scLMM/LMM-scRNAseq//Data/PBMC_Lupus_Kang8vs8_data_counts.rds")
dim(pbmc)
pbmc = pbmc[rowSums(pbmc)!=0,colSums(pbmc)!=0]
dim(pbmc)
Y_pbmc = GetAssayData(pbmc)
Y_pbmc[100:120,100:150]
nGenes <- colSums(Y_pbmc)
Ylog_pbmc = log2(t(Y_pbmc) + 1)
Ylog_pbmc <- as.matrix(Ylog_pbmc)
head(Ylog_pbmc[1:10,1:10])
head(Y[1:10,1:10])

##design matrix for random effects
pbmc$cell[is.na(pbmc$cell)] = 'unknown'
table(pbmc$cell)
sum(table(pbmc$cell)) == ncol(pbmc)
sum(table(pbmc$ind)) == ncol(pbmc)
Z_pbmc = data.frame(ind=pbmc$ind)
Z_pbmc = model.matrix(~0+ind, Z_pbmc) ### double check this with changjiang
#Z_pbmc = model.matrix(~0+ind, Z_pbmc) ### double check this with changjiang

head(Z_pbmc,3)

table(pbmc$cell)
##dimemsion of random effects
d_pbmc <- ncol(Z_pbmc)

##design matrix for fixed effects
X_pbmc = data.frame(stim=pbmc$stim,cell=pbmc$cell)
X_pbmc = model.matrix(~ 0 + cell + stim:cell, X_pbmc)

X_pbmc <- cbind(X_pbmc, log_nGenes = log(nGenes))
#colnames(X_pbmc)[1] = 'Intercept'
head(X_pbmc)
class(X_pbmc)


dim(X_pbmc)
dim(Z_pbmc)
dim(Ylog_pbmc)


##############################
## **Fit LMM by lmmfit.**
#########################
SEED <- 57774
set.seed(SEED)
s0 <- c(runif(length(d_pbmc)), 1)
maxIter <- 200
t1 <- Sys.time()
fit <- NULL
fit <- lmmfit(Y = Ylog_pbmc, X = X_pbmc, Z = Z_pbmc, 
              d = d_pbmc, max.iter = maxIter, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1) 

str(fit)
table(fit$niter)
sum(abs(fit$dlogL) > 1e-3)
dim(Ylog_pbmc)

##number of convergences
sum(fit$niter < maxIter) #18714
##number of non-convergences at epsilon = 1e-5
sum(fit$niter == maxIter) #176
##at epsilon = 1e-2
sum(apply(abs(fit$dlogL) > 1e-2, 2, any)) #153

##variance components of random effects
k <- 1 
range(fit$theta[k,]) # -4.487325e-05  2.136005e-01

##p-values for hypothesis tests of variance components: ##H0: theta <=0 vs H1: theta > 0
p <- pnorm(fit$theta[k, ]/fit$se[k, ], lower.tail = F) 
sum(p <= 0.05) # 854
range(p) #  0.0307845 1.0000000

##LMM tests
test <- lmmtest(fit)
test <- test[, grep("stim", colnames(test))]
test <- as.data.frame(test)
test$stim_FDR <- p.adjust(test$`cellCD8 T cells:stimstim_pvalue`, method = "BH") 
test$stim_FDR <- p.adjust(test$`cellB cells:stimstim_pvalue`, method = "BH") 

test <- test[order(test$stim_FDR),]
sum(test$stim_FDR <= 0.05, na.rm = T)

pv <- test[, grep("pvalue", colnames(test))]
par(mfrow = c(2,1), mar = c(4.5, 4.5, 1.1, 2.1))
qqpvalue(pv, col = "blue", cex = 0.6)
hist(pv$`cellCD8 T cells:stimstim_pvalue`, xlab = "Model I p-values for cellCD8 T cells:stim", col = "green", cex = 0.6, main=NA)


par(mfrow = c(2,5), mar = c(5.1, 4.1, 3.1, 1.1)) 
for (i in 1:ncol(pv)) {
  main <- gsub("_pvalue", "-wildtype", colnames(pv)[i])
  qqpvalue(pv[,i], col = "blue", cex = 0.6, main = main, cex.main = 0.8) }

dev.off()
par(mfrow = c(2,5), mar = c(5.1, 4.1, 3.1, 1.1)) 
for (i in 1:ncol(pv)) {
  main <- gsub("_pvalue", "-wildtype", colnames(pv)[i])
  hist(pv[,i], xlab = "p-values", col = "green", cex = 0.6, main = main, cex.main = 0.8) }

##############################
## Fit LMM by lmmfitSS.**
##############################

X = X_pbmc; Y= Ylog_pbmc; Z= Z_pbmc; d=d_pbmc
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

#save(fitss, lmm, t1, t2, file = paste0(dirOut, "/rat_set1_lmm.RData"))

pv <- lmm[, grep("pvalue", colnames(lmm))]
pv <- pv[, grep("stim", colnames(pv))]
head(pv)

hist(c(pv))

lmm_stim <- data.frame(lmm[, grep("stim", colnames(lmm))])
head(lmm_stim)
lmm_stim$score = -log(lmm_stim$stimstim_pvalue+1e-10)*lmm_stim$stimstim_t
lmm_stim = lmm_stim[order(lmm_stim$score, decreasing = T),]
lmm_stim = lmm_stim[!is.na(lmm_stim$stimstim_t),]
head(lmm_stim, 20)
tail(lmm_stim, 20)

head(rownames(lmm_stim), 20)
tail(lmm_stim, 20)


#####
##Test contrasts.
##Specify contrasts of interest.
nC <- 2
C <- matrix(0, nrow(fitss$coef), nC)
rownames(C) <- rownames(fitss$coef)
colnames(C) <- 1:nC
##cluster2 vs cluster1
j <- 1
C["cluster1", j] <- -1
C["cluster2", j] <- 1
colnames(C)[j] <- "cluster2-cluster1"

##cluster3 vs cluster1
j <- 2
C["cluster1", j] <- -1
C["cluster3", j] <- 1
colnames(C)[j] <- "cluster3-cluster1"
C

test <- lmmtest(fitss, contrast = C)
head(test)

