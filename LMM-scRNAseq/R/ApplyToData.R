source("~/scLMM/LMM-scRNAseq/lmmfit.R")
source("~/scLMM/LMM-scRNAseq/lmmfitSS.R")
source("~/scLMM/LMM-scRNAseq/lmmtest.R")
source("~/scLMM/LMM-scRNAseq/qqpvalue.R")

dirWork = '~/scLMM/LMM-scRNAseq/'
setwd(dirWork)
dir(dirWork)

library(Matrix)
library(MASS)
BiocManager::install('maditr')
library(maditr)
library(data.table)

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
Z_pbmc = data.frame(cell=pbmc$cell,ind=pbmc$ind)
Z_pbmc = model.matrix(~0+ cell+ind, Z_pbmc) ### double check this with changjiang
#Z_pbmc = model.matrix(~0+ind, Z_pbmc) ### double check this with changjiang

head(Z_pbmc,3)
head(Z)
table(pbmc$cell)
##dimemsion of random effects
d_pbmc <- ncol(Z_pbmc)

##design matrix for fixed effects
X_pbmc = data.frame(stim=pbmc$stim)
X_pbmc = model.matrix(~ stim, X_pbmc)
X_pbmc <- cbind(X_pbmc, log_nGenes = log(nGenes))
colnames(X_pbmc)[1] = 'Intercept'
head(X_pbmc)
head(X)
class(X)
class(X_pbmc)


dim(X_pbmc)
dim(Z_pbmc)
dim(Ylog_pbmc)

dim(X)
dim(Z)
dim(Y)

##############################
## **Fit LMM by lmmfit.**
#########################
SEED <- 57774
set.seed(SEED)
s0 <- c(runif(length(d_pbmc)), 1)

t1 <- Sys.time()
fit <- NULL
fit <- lmmfit(Y = Ylog_pbmc, X = X_pbmc, Z = Z_pbmc, 
              d = d_pbmc, s0 = s0, max.iter = 50, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1) 

str(fit)
table(fit$niter)
sum(abs(fit$dlogL) > 1e-3)
dim(Ylog_pbmc)


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
dim(lmm)
head(lmm)
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




##load data
##Y: counts
load("LMM-scRNAseq/rat_set1_XYZ.RData")
head(X)
dim(X)
head(Y) ### count
dim(Y)
head(Z)
dim(Z)
class(X)
class(Z)
class(Y)

##number of genes per cell
nGenes <- colSums(Y)
rm(X, Y, Z)

##Y: log2(t(counts) + 1)
load("LMM-scRNAseq/rat_set1_XZlogY.RData")
head(X)
head(Y)  ### log2(t(counts) + 1)
head(Z)

##design matrix for fixed effects
##cell clusters
dim(X)
X <- cbind(X, log_nGenes = log(nGenes))
head(X)

##design matrix for random effects
strainLew <- as.factor(Z)
Z <- model.matrix(~ 0 + strainLew)
head(Z)

##dimemsion of random effects
d <- ncol(Z)

##log2(t(counts) + 1)
dim(Y)
##Operating on "matrix" or "array"  is faster than "dgCMatrix"!!!
Y <- as.matrix(Y)
head(Y[1:10,1:10])



##############################
## **Fit LMM by lmmfit.**
#########################
SEED <- 57774
set.seed(SEED)
s0 <- c(runif(length(d)), 1)

t1 <- Sys.time()
fit <- NULL
fit <- lmmfit(Y = Y, X = X, Z = Z, d = d, s0 = s0, max.iter = 50, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1) 

str(fit)
table(fit$niter)
sum(abs(fit$dlogL) > 1e-3)

##############################
## Fit LMM by lmmfitSS.**
##############################

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
fitss <- lmmfitSS(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, Ynorm = Ynorm, n = n, d = d, 
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

#####
##Test all coefficients.
##cluster k vs cluster 10 (intercept)
lmm <- lmmtest(fitss)
dim(lmm)
head(lmm)
#save(fitss, lmm, t1, t2, file = paste0(dirOut, "/rat_set1_lmm.RData"))

pv <- lmm[, grep("pvalue", colnames(lmm))]
pv <- pv[, grep("cluster", colnames(pv))]
head(pv)

hist(c(pv))


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


##############################
### Comparison with NEBULA


##Load the NEBULA results.
load("rat_set1_nebula.RData")
difftime(t2, t1)

##nebula outputs
##summary (statistics): 
##The estimated coefficient, standard error and p-value for each predictor.
str(negbn)
st <- negbn$summary
rownames(st) <- st$gene
head(st)

##The genes which the LMM fit wasn't convergent.
table(fitss$niter)
indexNotconverge <- (fitss$niter >= 50)
sum(indexNotconverge)

##p-values
##NEBULA
pnb <- as.matrix(st[, grep("p_", colnames(st))])
head(pnb)

##LMM
plmm <- lmm[, grep("_p", colnames(lmm))]
plmm <- plmm[, -grep("log_", colnames(plmm))]
head(plmm)

i <- (!indexNotconverge)

par(mfrow = c(2, 1), mar = c(4.1,4.1,1.1,1.1))
hist(plmm[i, -1], xlab = "LMM p-values", main = NA)
hist(pnb[i, -1], xlab = "NEBULA pvalue", main = NA)

##QQ-plot
par(mfrow = c(2, 1), mar = c(4.1,4.1,1.1,1.1))
qqpvalue(plmm[i, -1], main = "LMM p-values", cex.lab = 0.8, cex.main = 0.8)
qqpvalue(pnb[i, -1], main = "NEBULA pvalue", cex.lab = 0.8, cex.main = 0.8)

##number of significant DE genes
##LMM
sum(plmm[i, -1] <= 0.05/sum(i), na.rm = T)

##NEBULA
sum(pnb[i, -1] <= 0.05/sum(i))


#####
##coefficients
bnb <- as.matrix(st[, grep("logFC", colnames(st))])
head(bnb)
dim(bnb)

blmm <- lmm[, grep("_t", colnames(lmm))]
blmm <- blmm[, -grep("log_", colnames(blmm))]
head(blmm)
dim(blmm)

i <- (!indexNotconverge)
par(mfrow = c(3,3), mar = c(4.1,4.1,1.1,1.1))
for (j in 2:ncol(blmm)) {
  plot(bnb[i, j], blmm[i, j], ylab = "LMM coef", xlab = "NEBULA coef")
  mtext(gsub("_.*", "", colnames(plmm)[j]), line = -1, cex = 0.8)
}


