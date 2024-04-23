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
#XY <- t(Y%*%X) #  argument is not a matrix
XY <- t(as.matrix(Y%*%X))
#ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y) #  argument is not a matrix
ZY <- t(as.matrix(Y%*%Z)) #as.matrix(t(Z)%*%Y) #  argument is not a matrix


ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)

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

print(paste0('time: ', difftime(t2, t1)))
# saveRDS(fitss, '~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021_X_sexAgeSeason_Z_donor.rds')
fit <- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021_X_sexAgeSeason_Z_donor.rds')
# Time difference of 1.200262 hours

## Time difference of 6.184403 hours
fit<- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021.rds')


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


colnames(plmm)
###### p value dataframe
a_cell_type_control = "CD4p_activated:CONTROL_pvalue"
df = data.frame(genes=rownames(plmm),cell_TB_pval=plmm[,a_cell_type_control])
head(df)
summary(df$cell_TB_pval)
df = df[order(df$cell_TB_pval, decreasing=F),]
head(df,30)

##### t value dataframe
df = data.frame(genes=rownames(plmm),cell_TB_tval=plmm[,a_cell_type_control])
head(df)
summary(df$cell_TB_tval)
df = df[order(df$cell_TB_tval, decreasing=T),]
head(df)

#####
colnames(test)
a_cell_type_control = 'CD4p_activated:CONTROL'#'Proximal_Tubule:Male'
a_cell_type_control_df = data.frame(genes= rownames(test),test[,grep(a_cell_type_control, colnames(test))])
head(a_cell_type_control_df)
a_cell_type_control = gsub(':', '.', a_cell_type_control)
a_cell_type_control_df$score = -log10(a_cell_type_control_df[[paste0(a_cell_type_control, '_pvalue')]])* abs(a_cell_type_control_df[[paste0(a_cell_type_control, '_t')]])
a_cell_type_control_df = a_cell_type_control_df[order(a_cell_type_control_df$score, decreasing = T),]
head(a_cell_type_control_df,30)
a_cell_type_control_df_control = a_cell_type_control_df[a_cell_type_control_df[[paste0(a_cell_type_control, '_t')]]>0,]
head(a_cell_type_control_df_control, 20)
dev.off()
gridExtra::grid.table(head(a_cell_type_control_df_control, 20))

a_cell_type_control_df_TB = a_cell_type_control_df[a_cell_type_control_df[[paste0(a_cell_type_control, '_t')]]<0,]
head(a_cell_type_control_df_TB, 20)
dev.off()
gridExtra::grid.table(head(a_cell_type_control_df_TB, 20))
getwd()



source('~/RatLiver/Codes/Functions.R')
Initialize()
library(gprofiler2)
library(ggplot2)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
head(a_cell_type_control_df_male,30)
num_genes = 200

enrich_res = get_gprofiler_enrich(markers=a_cell_type_control_df_control$genes[1:num_genes], model_animal_name)
head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
title = a_cell_type_control
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))


num_genes = 200
enrich_res = get_gprofiler_enrich(markers=a_cell_type_control_df_TB$genes[1:num_genes], model_animal_name)
head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
enrich_res_pos = enrich_res_pos[!is.na(enrich_res_pos$log_p),]
title = gsub(pattern = 'CONTROL', 'TB', a_cell_type_control)
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))


table(fit$niter)



#####
library(grid); 
start_i = 101#93
colnames_to_check = colnames(test)[start_i:length(colnames(test))]
a_cell_type_control = gsub('_pvalue','',colnames_to_check)[i]
i = 3
num_genes = 300
#pdf('~/scLMM/LMM-scRNAseq-jan2024/result_tables_DEs_Nathan_NatImm_2021_V2.pdf', width = 14, height = 9)
pdf('~/scLMM/LMM-scRNAseq-jan2024/result_tables_DEs_Nathan_NatImm_X_sexAgeSeason_Z_donor.pdf', width = 14, height = 9)
for(i in 1:length(colnames_to_check)){ #
  a_cell_type_control = gsub('_pvalue','',colnames_to_check)[i]
  a_cell_type_control_df = data.frame(genes= rownames(test),test[,grep(a_cell_type_control, colnames(test))])
  head(a_cell_type_control_df)
  a_cell_type_control = gsub(':', '.', a_cell_type_control)
  a_cell_type_control = gsub(':', '.', a_cell_type_control)
  #a_cell_type_control_df$score = -log10(a_cell_type_control_df[[paste0(a_cell_type_control, '_pvalue')]])* abs(a_cell_type_control_df[[paste0(a_cell_type_control, '_t')]])
  a_cell_type_control_df$score = -log10(a_cell_type_control_df[[3]])* abs(a_cell_type_control_df[[2]])
  a_cell_type_control_df = a_cell_type_control_df[order(a_cell_type_control_df$score, decreasing = T),]
  head(a_cell_type_control_df,30)
  a_cell_type_control_df_control = a_cell_type_control_df[a_cell_type_control_df[[2]]>0,]
  head(a_cell_type_control_df_control, 30)
  #dev.off()
  grid.newpage()
  gridExtra::grid.table(head(a_cell_type_control_df_control, 20))
  
  
  enrich_res = get_gprofiler_enrich(markers=a_cell_type_control_df_control$genes[1:num_genes], model_animal_name)
  if(!is.null(enrich_res)){
    print('YES!')
    head(enrich_res$result,30)
    enrich_res_pos = data.frame(enrich_res$result)
    enrich_res_pos = enrich_res_pos[1:20,]
    enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
    enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
    title = a_cell_type_control
    p=ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
      theme_classic()+ylab('')+ggtitle(paste0(title))
    print(p)
    
  }
  
  
  a_cell_type_control_df_TB = a_cell_type_control_df[a_cell_type_control_df[[2]]<0,]
  head(a_cell_type_control_df_TB, 20)
  #dev.off()
  grid.newpage()
  gridExtra::grid.table(head(a_cell_type_control_df_TB, 20))
  
  
  enrich_res = get_gprofiler_enrich(markers=a_cell_type_control_df_TB$genes[1:num_genes], model_animal_name)
  if(!is.null(enrich_res)){
    print('YES!')
    head(enrich_res$result,30)
    enrich_res_pos = data.frame(enrich_res$result)
    enrich_res_pos = enrich_res_pos[1:20,]
    enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
    enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
    enrich_res_pos = enrich_res_pos[!is.na(enrich_res_pos$log_p),]
    title = gsub(pattern = 'CONTROL', 'TB', a_cell_type_control)
    p=ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
      theme_classic()+ylab('')+ggtitle(paste0(title))
    print(p)
    
  }
  
  
}
dev.off()




##histograms of p-values


tvlmm <- test[, grep("_t", colnames(test)), drop = F]
i <- (apply(!is.na(tvlmm), 1, all)) 
jset <- grep(":CONTROL", colnames(tvlmm))
length(jset)
sum(i)
dim(tvlmm)
length(jset)
pdf('~/scLMM/LMM-scRNAseq-jan2024/results_pvalues_Nathan_NatImm_2021.pdf')
pdf('~/scLMM/LMM-scRNAseq-jan2024/results_pvalues_Nathan_NatImm_X_sexAgeSeason_Z_dono.pdf')

par(mfrow = c(4, 4))
for (j in jset){
  #h1 <- hist(pv[i,j], plot = F)
  #h2 <- hist(plmm[i,j], plot = F)
  #ylim <- c(0, 1.1*max(h2$counts, h2$counts))
  nm <- gsub("p", "+", gsub("_pvalue", "", colnames(plmm)[j]))
  
  #par(mar = c(5.1,4.1,2.1,0))
  #nm <- gsub("p_", "", colnames(pv)[j])
  #hist(pv[i,j], ylim = ylim, xlab = "nebula p-values", 
  #     main = nm, cex.main = 0.8, cex.axis = 0.8)
  #par(mar = c(5.1,2.1,2.1,2.1))	
  #nm <- gsub("_pvalue", "", colnames(plmm)[j])
  hist(plmm[i, j], xlab = "lmmfit p-value", ylab = NULL, 
       #ylim =ylim,
       main = nm, cex.main = 0.8, cex.axis = 0.8, yaxt = "n")
}

dev.off()
