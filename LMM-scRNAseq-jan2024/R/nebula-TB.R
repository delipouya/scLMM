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
table(coldata$age)
table(coldata$season)

summary(coldata$prop_NAT)


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

##raw counts
Y <- counts
dim(Y) 
#[1]  7017 26820

nGenes <- colSums(Y)
##
rm(counts)
table(colSums(table(coldata$batch, coldata$donor) >0))
colnames(coldata)
gc()



##nebula
##fixed effect desigm matrix

coldata$cluster_name[is.na(coldata$cluster_name)] = 'unknown'
X <- model.matrix(~ log(nGene) + sex + age + season + cluster_name + cluster_name:TB_status, data = coldata)
colnames(X)
colnames(X) <- gsub("cluster_name", "", colnames(X))
colnames(X) <- gsub("TB_status", "", colnames(X))
#colnames(X) <- gsub("sex", "", colnames(X))
#colnames(X) <- gsub("age", "", colnames(X))
#colnames(X) <- gsub("season", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

head(X)
dim(X)


##random effect (sample groups)
Z <- data.frame(donor=coldata$donor)
Z <- model.matrix(~ 0 +donor,data = coldata)


##########################################################################################

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



##random effect (sample groups)
Z <- coldata$donor
table(Z)
length(Z) #[1] 26820


##model:
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model (the same model as that in the lme4 package).

model <- "NBLMM" #model = "NBGMM",

t1 <- Sys.time()
negbn <- NULL
##The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = Y[, o], id = as.factor(Z[o]), pred = X[o, ], 
                model = model, covariance = T, output_re = T)
t2 <- Sys.time()

difftime(t2, t1)
#Time difference of 53.49769 mins

timeUnits <- "secs"
rtnebula <- difftime(t2, t1, units = timeUnits)
print(paste0('time difference(min) is: ', difftime(t2, t1))) 
#"time difference(min) is: 28.1434102137883"
print(paste0('time difference(sec) is: ', rtnebula)) 
# "time difference(sec) is: 1688.6046128273"

##nebula outputs
##summary (statistics): 
##The estimated coefficient, standard error and p-value for each predictor.
str(negbn)

saveRDS(negbn, 'nebula_TB_allCovs_interaction.rds')
#negbn = readRDS('nebula_stimPBMC.rds')



plot(negbn$overdispersion$Subject)
plot(negbn$overdispersion$Cell)

##
##https://github.com/lhe17/nebula
##  1: The convergence is reached due to a sufficiently small improvement of the function value.
##-10: The convergence is reached because the gradients are close to zero 
##     (i.e., the critical point) and no improvement of the function value can be found.
table(negbn$convergence)
# -40  -10    1 
#   5  101 6911 


##################################################