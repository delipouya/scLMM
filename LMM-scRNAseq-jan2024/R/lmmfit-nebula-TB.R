##########
##Fit LMM by lmmfit
#source("~/scLMM/LMM-scRNAseq/R/lmmfit..R")
source("~/scLMM/LMM-scRNAseq-jan2024//R/lmmfit.nt.R")
source("~/scLMM/LMM-scRNAseq-jan2024//R/lmmfitSS.R")
source("~/scLMM/LMM-scRNAseq-jan2024//R/lmmtest.R")
source("~/scLMM/LMM-scRNAseq-jan2024//R/qqpvalue.R")

library(nebula)
library(Seurat)
library(MASS)

counts = readRDS('~/sciFA/Data/Nathan_NatImm_2021.rds')
coldata = counts@meta.data
summary(coldata)
colnames(coldata)
length(table(coldata$batch))
length(table(coldata$donor))
table(coldata$sex)


##data
dim(counts)
#[1]  33569 500089
counts = GetAssayData(counts)
##Filtering genes
##Mixed models in muscat package - 3.4 Cell-level analysis:
##(1) subpopulations with at least 10 cells in at least 2 samples (not necessary?)
##(2) genes with a count >= 1 in at least 20 cells

all(colnames(data) == rownames(coldata))
##number of celss
nCells <- rowSums(counts > 0)
hist(log2(nCells))

minCells <- 2^5

##number of cells in a group_id
nCellsgrp <- do.call(cbind, 
                     tapply(1:ncol(counts), as.factor(coldata$donor), 
                            function(j) rowSums(counts[, j, drop = F] > 0))
)

minCellsgrp <- 10

##number of counts
ncts <- rowSums(counts)
#hist(log2(ncts)) #outliers in the upper tail

maxCounts <- 2^20
sum(ncts > maxCounts)

##number of counts in a group_id
nctsgrp <- do.call(cbind, 
                   tapply(1:ncol(counts), as.factor(coldata$donor), 
                          function(j) rowSums(counts[, j, drop = F]))
)


##nebula filtering:
##Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.

cpc <- rowSums(counts)/ncol(counts)
sum(cpc <= 0.005) 
#[1] 31


##Filtering
minCells <- 32 #2^5
minCellsgrp <- 15
maxCounts <- 2^20
minCountsgrp <- 2*minCellsgrp
minCountsgrp
mincpc <- 0.005

index <- (nCells >= minCells) & (rowSums(nCellsgrp >= minCellsgrp) >= 2)
index <- index & (ncts <= maxCounts) & (rowSums(nctsgrp >= minCountsgrp) >= 2)
sum(index)
index <- index & (cpc > mincpc)
sum(index)

counts <- counts[index, ] ### around 100 genes are removed
rm(index)
dim(counts)
#[1]  11322 500089


### >>> start here
##raw counts
Y <- counts
dim(Y) 
#[1]  7017 26820

nGenes <- colSums(Y)
##
rm(counts)
table(colSums(table(coldata$batch, coldata$donor) >0))
colnames(coldata)
##################################################
#NEBULA
##https://github.com/lhe17/nebula
##Checking convergence for the summary statistics and quality control
##  1: The convergence is reached due to a sufficiently small improvement of the function value.
##-10: The convergence is reached because the gradients are close to zero 
##     (i.e., the critical point) and no improvement of the function value can be found.
##
##Depending on the concrete application, 
##the estimated gene-specific overdispersions can also be taken into consideration in quality control. 
##For example, when testing differential expression for a variable, 
##genes with a very large estimated cell-level overdispersion should be filtered out because such genes have huge unexplained noises.
##
##If the variable of interest is subject-level, 
##genes with a very large subject-level overdispersion (>1) 
##should be removed or interpreted cautiously as well.
##
##The NBLMM is the same model as that adopted in the glmer.nb function in the lme4 R package, 
##but is computationally much more efficient by setting method='LN'. 


##nebula
##fixed effect desigm matrix
X <- model.matrix(~ log(nGene) + cluster_name + cluster_name:TB_status, data = coldata)
colnames(X)
colnames(X) <- gsub("cluster_name", "", colnames(X))
colnames(X) <- gsub("TB_status", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

head(X)
dim(X)

##random effect (sample groups)
Z <- data.frame(donor=coldata$donor, batch=coldata$batch)
Z <- model.matrix(~donor + batch,data = coldata)
table(Z)
length(Z) #[1] 26820
head(Z)



##################################################
##lmmfit
dim(Y)
#[1]  7017 26820

##log-transformation
##log2(1 + counts)
##log2(1+Y) by groupping to reduce data size
ngrp <- 3
sizegrp <- round(nrow(Y)/ngrp)
for (i in 1:ngrp){
  j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
  print(range(j))
  Y[j, ] <- log2(1 + Y[j, ])
}

##transpose
#Y <- t(Y)



dim(Y) 
#[1] 26820  7017

##fixed effect desigm matrix >> START
coldata$cluster_name[is.na(coldata$cluster_name)] = 'unknown'
X <- model.matrix(~ log(nGene) + cluster_name + cluster_name:TB_status, data = coldata)
colnames(X)
colnames(X) <- gsub("cluster_name", "", colnames(X))
colnames(X) <- gsub("TB_status", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))
head(X)
dim(X)

dim(coldata)
sum(table(coldata$cluster_name))
sum(table(coldata$TB_status))


##random effect (sample groups)
Z <- data.frame(donor=coldata$donor, batch=coldata$batch)
Z <- model.matrix(~ 0 +donor + batch,data = coldata)
colnames(Z)
d <- ncol(Z)
##########
##Fit LMM by lmmfit

##Operating on "matrix" or "array"  is faster than "dgCMatrix"!!!
#Y <- as.matrix(Y)

dim(Y)
dim(Z)
dim(X)
timeUnits <- "secs"
maxIter <- 200 
epsilon <- 1e-8 

################# ################# ################# 
########################## lmmfit.SS
XY <- t(Y%*%X) #  argument is not a matrix
#XY <- t(as.matrix(Y%*%X))
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y) #  argument is not a matrix
#ZY <- t(as.matrix(Y%*%Z)) #as.matrix(t(Z)%*%Y) #  argument is not a matrix
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

XXinv <- ginv(t(X)%*%X)
Ynorm <- rowSums(Y*Y) #colSums(Y*Y)
XY <- as.matrix(XY)
ZY <- as.matrix(ZY)
n = ncol(Y) ## sample size 
t1 <- Sys.time()
fitss <- lmmfitSS(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, 
                  Ynorm = Ynorm, n = n, d = d, max.iter = 100, epsilon = 1e-5)
t2 <- Sys.time()
difftime(t2, t1)



################# ################# ################# ################# 
################# performing the fit using lmmfit.nt 
t1 <- Sys.time()
fit <- NULL
fit <- lmmfit.nt(Y = Y, X = X, Z = Z, d = d, max.iter = maxIter, epsilon = epsilon)
t2 <- Sys.time()
rtlmm <- difftime(t2, t1, units = timeUnits) 

table(fit$niter)


##################################################
##comparison of results

##########
##lmmfit

rtlmm; timeUnits; maxIter; epsilon
#Time difference of 49.96773 secs
#[1] "secs"
#[1] 200
#[1] 1e-08

##fixed effects
felmm <- fit$coef

##variance components
slmm <- fit$theta

##LMM tests
test <- lmmtest(fit)
dim(test)
head(test)

##t-values
tvlmm <- test[, grep("_t", colnames(test)), drop = F]
dim(tvlmm)
head(tvlmm)

##p-values
plmm <- test[, grep("_pvalue", colnames(test)), drop = F]
dim(plmm)
head(plmm)