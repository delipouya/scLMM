source('~/RatLiver/Codes/Functions.R')
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

############
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
rm(your_scRNAseq_data_object)
gc()

sample_size = 3000

# add_info = '' 
# add_info = '_sample_RE_strain_RE'
# add_info = '_sample_RE_strain_FE'
# add_info = '_strain_RE_BeforeImpute'
add_info = '_strain_RE_AfterImpute'

#### making sure the 3000 subsample is a part of the 15000 sample set 
# umi_samples_15000 = readRDS('ratLiver_LMM_results/subsample_results/sample_15000/Rat_subsampleUMIs_samleSize_150002.rds')
# selected_umi_index = sample(size = sample_size, umi_samples_15000)

umi_samples_3000 = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/selected_umi_3000_outof_15000.rds')
merged_samples_sub = merged_samples[,umi_samples_3000]
#### Making the input data even smaller
merged_samples_sub = merged_samples_sub[1:500,1:1000]

################### impute the single cell data ################
######## replace the zeros with NA values
data.mis = as.matrix(GetAssayData(merged_samples_sub))
data.mis = data.mis[rowSums(data.mis) != 0, colSums(data.mis) != 0]
dim(data.mis)

hist(rowSums(data.mis), breaks = 100)
hist(colSums(data.mis))
data.mis[data.mis==0] <- NA

##### evaluating the number of NAs in the dataset prior to the imputation
num_NA_before_impute = unlist(lapply(data.frame(t(data.mis)), function(x) sum(is.na(x)) ))
hist(num_NA_before_impute,breaks = 500,
     main = 'Distribution of NA values\nfor each genes before imputation\ndataset: 500 genes, 1000 UMIs',
     xlab = 'Number of NAs in gene expression before imputation')

##### impute using the MICE package
imputed_Data <- mice(t(data.mis)) ## , diagnostics = FALSE , remove_collinear = FALSE
imputed_Data_2 <- mice(t(data.mis), remove_collinear = FALSE)  ### leads to error

completedData.withNA <- complete(imputed_Data,1)
dim(completedData.withNA)
head(completedData.withNA)
#saveRDS(completedData.withNA, '~/scLMM/ratLiver_LMM_results/Imputed_subset_RatData.rds') 
#saveRDS(t(data.mis), '~/scLMM/ratLiver_LMM_results/Before_impute_subset_RatData.rds') 

num_NA_imputed_genes = unlist(lapply(completedData.withNA, function(x) sum(is.na(x)) ))
hist(num_NA_imputed_genes,breaks = 500,
     main = 'Distribution of NA values\nfor each genes after imputation\ndataset: 500 genes, 1000 UMIs',
     xlab = 'Number of NAs in the imputated gene expression')

##### impute using the Hmisc package
library(Hmisc)
hmisc_imputed_data = Hmisc::impute(t(data.mis), 'random')
dim(hmisc_imputed_data)
lapply(hmisc_imputed_data, head)
hmisc_imputed_data.df = as.matrix(hmisc_imputed_data) #### how to work with the output?


##### impute using the missForest package
library("missForest")
missFor_imputed_data <- missForest(t(data.mis)) ### leads to an error

#library("eimpute")
#rank = qr(data.mis)$rank
#data_impute <- eimpute(data.mis, rank)
################################################################

############### reading the subsetted data before and after imputation #################
data_before_impute = readRDS( '~/scLMM/ratLiver_LMM_results/Before_impute_subset_RatData.rds') 
data_after_impute = readRDS( '~/scLMM/ratLiver_LMM_results/Imputed_subset_RatData.rds') 

data = as.matrix(t(data_after_impute))
data[is.na(data)] <- 0
dim(data)


