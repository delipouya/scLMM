
##################################################
##load data
##PBMC data
library(nebula)
library(Seurat)

dirData = '~/scLMM/LMM-scRNAseq/Data/'
datafile <- paste0(dirData, "/Human_Kidney_data.rds")
data = readRDS(file = datafile)

coldata = data@meta.data
counts = GetAssayData(data, assay = 'RNA', layer = 'counts')
rowSums(counts)

##'sampleID': unique sample identifiers
##'Cell_Types_Broad': subpopulation (cell cluster) assignments 
##'sex': make vs female experimental group/condition (control/treatment, healthy/diseased)

table(coldata$sampleID)
# CD45_1 CD45_10  CD45_2  CD45_3  CD45_4  CD45_5  CD45_6  CD45_7  CD45_8  CD45_9  Total1  Total2  Total3  Total4  Total5  Total6  Total7 
# 584     714     315     938     126     344    1102    1182     553    1041    1080    4738    2366    3739    3280     606    1811 
# Total8  Total9 
# 1894    1264 
table(coldata$sex)
# Female   Male 
# 15945  11732 
table(coldata$sampletype)
table(coldata$Cell_Types_Broad)
table(coldata$CMVstatus)


##counts
dim(counts)
#[1]  22484 27677

##Filtering genes
##Mixed models in muscat package - 3.4 Cell-level analysis:
##(1) subpopulations with at least 10 cells in at least 2 samples (not necessary?)
##(2) genes with a count >= 1 in at least 20 cells

all(colnames(counts) == rownames(coldata))

##number of celss
nCells <- rowSums(counts > 0)
hist(log2(nCells))

minCells <- 2^6

##number of cells in a group_id (sex)
nCellsgrp <- do.call(cbind, 
                     tapply(1:ncol(counts), as.factor(coldata$sex), 
                            function(j) rowSums(counts[, j, drop = F] > 0))
)

minCellsgrp <- 50#used to be 10

##number of counts
ncts <- rowSums(counts)
hist(log2(ncts)) #outliers in the upper tail

maxCounts <- 2^20
sum(ncts > maxCounts)

##number of counts in a group_id
nctsgrp <- do.call(cbind, 
                   tapply(1:ncol(counts), as.factor(coldata$sex), 
                          function(j) rowSums(counts[, j, drop = F]))
)


##nebula filtering:
##Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.

cpc <- rowSums(counts)/ncol(counts)
sum(cpc <= 0.005) 
#[1] 8585 ## a bit too high!


##Filtering
minCells <- 64#32 #2^5
minCellsgrp <- 50
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
#[1]  13884 27677 # minCellsgrp = 10
# 12628 27677



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


##raw counts
Y <- counts
dim(Y) 
#[1]  12628 27677

nGenes <- colSums(Y)

##
#rm(counts)

##nebula
##fixed effect desigm matrix
X <- model.matrix(~ log(nGenes) + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
colnames(X) <- gsub("Cell_Types_Broad", "", colnames(X))
colnames(X) <- gsub("sex", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

head(X)
dim(X)

##random effect (sample groups)
Z <- coldata$sampleID
table(Z)
length(Z) #[1] 27677


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
#Error in { : task 2138 failed - "objective in x0 returns NA"
#In addition: There were 50 or more warnings (use warnings() to see the first 50)
difftime(t2, t1)
#Time difference of 53.49769 mins

timeUnits <- "secs"
rtnebula <- difftime(t2, t1, units = timeUnits)
print(paste0('time difference(min) is: ', difftime(t2, t1))) # "time difference(min) is: 2.84918628692627"
print(paste0('time difference(sec) is: ', rtnebula)) #  "time difference(sec) is: 10257.0706329346"

#saveRDS(negbn, 'nebula_kidneyhuman.rds')
negbn = readRDS('~/scLMM/LMM-scRNAseq/R/nebula_kidneyhuman.rds')


##nebula outputs
##summary (statistics): 
##The estimated coefficient, standard error and p-value for each predictor.
str(negbn)

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
##lmmfit
dim(Y)
#[1]  7017 26820
## 12628 27677

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
Y <- t(Y)

dim(Y) 
#[1] 26820  7017

##random effects
##sample groups
Z <- model.matrix(~ 0 + sampleID, data = coldata)
colnames(Z) <- gsub(".*sampleID", "", colnames(Z))
dim(Z) #[1] 26820    8
head(Z)

d <- ncol(Z)


##########
##Fit LMM by lmmfit
source("~/scLMM/LMM-scRNAseq/R/lmmfit.R")
source("~/scLMM/LMM-scRNAseq/R/lmmfitSS.R")
source("~/scLMM/LMM-scRNAseq/R/lmmtest.R")
source("~/scLMM/LMM-scRNAseq/R/qqpvalue.R")

##Operating on "matrix" or "array"  is faster than "dgCMatrix"!!!
Y <- as.matrix(Y)

timeUnits <- "secs"
maxIter <- 200 
epsilon <- 1e-8 

t1 <- Sys.time()
fit <- NULL
fit <- lmmfit(Y = Y, X = X, Z = Z, d = d, max.iter = maxIter, epsilon = epsilon)
t2 <- Sys.time()
rtlmm <- difftime(t2, t1, units = timeUnits) 
# Time difference of 109.1035 secs
# Time difference of 91.97659 secs
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


##########
##nebula
rtnebula = NA
rtnebula[[1]] = 10257.0706329346
#Time difference of 3209.861 secs
rtnebula[[1]]/rtlmm[[1]]
#[1] 64.23869
str(negbn)

table(negbn$convergence)
# -40  -10    1 
#   5  101 6911 

#  -60   -40   -25   -10     1 
#11685    57   298     3   585 

st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

any(is.na(st))
sum(is.na(st))
st[rowSums(is.na(st)) > 0, ]
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

##fixed effects, se, and p-values
iFC <- grep("logFC_", colnames(st))
ise <- grep("se_", colnames(st))
ipv <- setdiff(grep("p_", colnames(st)), c(iFC, ise))

##fixed effects
b <- st[, iFC]

##se
se <- st[, ise]

##t-values
tv <- b/se

##p-values
pv <- st[, ipv]

range(pv - 2*pnorm(-abs(tv)), na.rm = T)
#[1] -2.331468e-15  2.386980e-15


##########
##comparison of lmmfit and nebula
##for the genes that lmmfit and nebula are convergent.

all(rownames(tv) == rownames(tvlmm))
any(is.na(tv))
any(is.na(tvlmm))


hist(negbn$overdispersion$Subject)
hist(negbn$overdispersion$Cell)
##genes with convergence
i <- (negbn$overdispersion$Subject < 4) ## used to be 1 
sum(i)
i <- i & (negbn$overdispersion$Cell < 100) & (negbn$convergence >= -10)
sum(i)
#[1] 7000
i <- i & (colSums(abs(fit$dlogL) < 1e-5) == nrow(fit$dlogL))
sum(i)
#[1] 6967
i <- i & (apply(!is.na(tv), 1, all)) & (apply(!is.na(tvlmm), 1, all))
sum(i)
#[1] 6958


##t-values
j <- 2:ncol(tv)
plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]), 
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1)

##p-values
j <- 2:ncol(pv)
#j <- grep(":stim", colnames(pv))
plot(as.matrix(-log10(plmm[i, j])), as.matrix(-log10(pv[i, j])), 
     xlab = "lmmfit -log10(p-values)", ylab = "nebula -log10(p-values)", cex = 0.6)
abline(0, 1)


old.par <- par(no.readonly = TRUE)
on.exit(par(old.par))
##cell-type specific
jset <- grep(":Male", colnames(pv))
nc <- round(sqrt(length(jset)))
nr <- ceiling(length(jset)/nc)
par(mfrow = c(nr, nc), mar = c(5.1,4.1,2.1,1.1))
for (j in jset){
  nm <- gsub("p", "+", gsub("_t", "", colnames(tvlmm)[j]))
  plot(tvlmm[i,j], tv[i,j], cex = 0.6,
       xlab = "lmmfit t-values", ylab = "nebula t-values", main = nm, cex.main = 0.8)
  abline(0,1)
}


##histograms of p-values
jset <- grep(":Male", colnames(pv))
length(jset)
par(mfrow = c(4, 4))
for (j in jset){
  h1 <- hist(pv[i,j], plot = F)
  h2 <- hist(plmm[i,j], plot = F)
  ylim <- c(0, 1.1*max(h1$counts, h2$counts))
  nm <- gsub("p", "+", gsub("_pvalue", "", colnames(plmm)[j]))
  
  par(mar = c(5.1,4.1,2.1,0))
  #nm <- gsub("p_", "", colnames(pv)[j])
  hist(pv[i,j], ylim = ylim, xlab = "nebula p-values", 
       main = nm, cex.main = 0.8, cex.axis = 0.8)
  
  par(mar = c(5.1,2.1,2.1,2.1))	
  #nm <- gsub("_pvalue", "", colnames(plmm)[j])
  hist(plmm[i, j], ylim = ylim, xlab = "lmmfit p-value", ylab = NULL, 
       main = nm, cex.main = 0.8, cex.axis = 0.8, yaxt = "n")
}



##barplot of common genes in the top
all(rownames(plmm) == rownames(pv))
glist <- rownames(plmm)[i]

jset <- grep(":Male", colnames(pv))
ntops <- c(250, 500)

commonGenes <- NULL
for (nTop in ntops){
  cgenes <- NULL	
  #par(mfrow = c(4, 4))
  for (j in jset){
    p <- NULL
    p <- pv[i,j]
    g1 <- glist[order(p)[1:nTop]]
    
    p <- NULL
    p <- plmm[i,j]
    g2 <- glist[order(p)[1:nTop]]
    
    ##proportion of the common genes in the top
    cgenes <- c(cgenes, length(intersect(g1, g2))/nTop)
  }
  commonGenes <- cbind(commonGenes, cgenes)
}
rownames(commonGenes) <- gsub("p_", "+_", gsub(":Male_pvalue", "", colnames(plmm)[jset]))
colnames(commonGenes) <- ntops
commonGenes
dev.off()
gridExtra::grid.table(commonGenes)

colset <- rainbow(nrow(commonGenes))

par(mar = c(5.1, 4.1, 4.1, 8.1))
barplot(commonGenes, beside = T, ylim = c(0, 1), col = colset,
        xlab = "Number of top genes", ylab = "Proportion of the common genes in the top")
legend("topleft", rownames(commonGenes), inset = c(1, 0), xpd = T, pch = 15, 
       col = colset, bty = "n", cex = 0.7)

##################################################

##### interpretation of the results

##fixed effects
felmm <- fit$coef
dim(felmm)
felmm = data.frame(t(felmm))
head(felmm)
##variance components
slmm <- fit$theta
slmm = data.frame(t(slmm))
head(slmm)
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
colnames(plmm)
###### p value dataframe
a_cell_type_sex = "Proximal_Tubule:Male_pvalue"
df = data.frame(genes=rownames(plmm),cell_male_pval=plmm[,a_cell_type_sex])
head(df)
summary(df$cell_male_pval)
df = df[order(df$cell_male_pval, decreasing=F),]
head(df,30)

##### t value dataframe
df = data.frame(genes=rownames(plmm),cell_male_pval=plmm[,a_cell_type_sex])
head(df)
summary(df$cell_male_pval)
df = df[order(df$cell_male_pval, decreasing=F),]

#####
a_cell_type_sex ='Proximal_Tubule:Male'
a_cell_type_sex_df = data.frame(genes= rownames(test),test[,grep(a_cell_type_sex, colnames(test))])
head(a_cell_type_sex_df)
a_cell_type_sex_df$score = -log10(a_cell_type_sex_df$Proximal_Tubule.Male_pvalue)* abs(a_cell_type_sex_df$Proximal_Tubule.Male_t)
a_cell_type_sex_df = a_cell_type_sex_df[order(a_cell_type_sex_df$score, decreasing = T),]
head(a_cell_type_sex_df,30)
a_cell_type_sex_df_male = a_cell_type_sex_df[a_cell_type_sex_df$Proximal_Tubule.Male_t>0,]
head(a_cell_type_sex_df_male, 16)
a_cell_type_sex_df_female = a_cell_type_sex_df[a_cell_type_sex_df$Proximal_Tubule.Male_t<0,]
head(a_cell_type_sex_df_female, 16)
