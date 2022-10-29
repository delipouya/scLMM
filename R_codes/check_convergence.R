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

library(ggvenn)
####### Functions

get_total_cells_expressed <- function(genes, UMIs ){
  rowSums(GetAssayData(merged_samples[genes, UMIs])!=0)
}
get_total_cells_empty <- function(genes, UMIs ){
  rowSums(GetAssayData(merged_samples[genes, UMIs])==0)
}
get_mean_expression <- function(genes, UMIs){
  rowMeans(GetAssayData(merged_samples[genes, UMIs]))
}
get_variance <- function(genes, UMIs){
  rowVars(as.matrix(GetAssayData(merged_samples[genes, UMIs])))
}

get_convergence_df <- function(genes, UMIs, tag=NULL){
  #### info on the converged genes  
  gene_info_df = data.frame(total_cells_expressed = get_total_cells_expressed(genes, UMIs), 
                            total_cells_empty = get_total_cells_empty(genes, UMIs),
                            mean_expression = get_mean_expression(genes, UMIs), 
                            variance = get_variance(genes, UMIs))
  
  if(length(tag)>0) {
    gene_info_df$tag = tag
  }
  
  return(gene_info_df)
}


get_plots <- function(conv_df, run_info, is_model_converged, perc_non_conv){
  
  p1=ggplot(conv_df, aes(x=total_cells_expressed, fill=tag)) +
    geom_histogram( color="black", alpha = 0.5, bins = 50, position = 'identity') +theme_classic() +
    labs(title = 'Total number of cells expressing the gene', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='total cells expressed', y='Count')
  
  p2=ggplot(conv_df, aes(x=tag, y=total_cells_expressed, fill=tag)) +
    geom_boxplot( color="black", alpha = 0.8, position = 'identity') +theme_classic() +
    labs(title = 'Total number of cells expressing the gene', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='',y='total cells expressed')
  
  p3=ggplot(conv_df, aes(x=total_cells_empty, fill=tag)) +
    geom_histogram( color="black", alpha = 0.5, bins = 50, position = 'identity') +theme_classic() +
    labs(title = 'Total number of empty cells', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='total number of empty cells', y='Count')
  
  p4=ggplot(conv_df, aes(x=tag, y=total_cells_empty, fill=tag)) +
    geom_boxplot( color="black", alpha = 0.8, position = 'identity') +theme_classic() +
    labs(title = 'Total number of empty cells', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='',y='total empty cells')
  
  p5=ggplot(conv_df, aes(x=tag,y=mean_expression, fill=tag)) +
    geom_boxplot( color="black", alpha = 0.5, position = 'identity') +theme_classic() +
    labs(title = 'The mean gene expression distribution', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='',y='mean gene expression')
  
  p6=ggplot(conv_df, aes(x=tag,y=variance, fill=tag)) +
    geom_boxplot( color="black", alpha = 0.5, position = 'identity') +theme_classic() +
    labs(title = 'The gene expression variance distribution', 
         subtitle = paste0(run_info,'\n' ,
                           sum(!is_model_converged), ' of total ', length(is_model_converged) , 
                           ' genes (=',perc_non_conv,'%) did not converge'), 
         x='',y='variance of gene expression')
  
  return(list(p1,p2,p3,p4,p5,p6))
  
}

#########  Import the data
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)


######### 15000 - strain as random effect
is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_15000/Rat_is_model_converged_samleSize_15000.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_1500 = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

UMIs_15000 = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_15000/Rat_subsampleUMIs_samleSize_150002.rds') 
converged_genes_15000 = names(is_model_converged)[is_model_converged]
Nonconverged_genes_15000 = names(is_model_converged)[!is_model_converged]

perc_non_conv_1500
length(Nonconverged_genes_15000)/(length(Nonconverged_genes_15000)+length(converged_genes_15000))


conv_df_15000 = rbind(get_convergence_df(converged_genes_15000, UMIs_15000, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_15000, UMIs_15000, tag='non_conv'))
head(conv_df_15000)
sum(conv_df_15000$total_cells_expressed <60)

plots_15000 = get_plots(conv_df_15000, run_info = '15000 cells - strain as random effect', 
                        is_model_converged =is_model_converged, 
                        perc_non_conv=perc_non_conv_1500)
plots_15000[[1]]
rm(is_model_converged)
############# 3000 strain as random effect 

is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000/Rat_is_model_converged_samleSize_3000_UMIs.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_3000 = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

UMIs_3000 = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/selected_umi_3000_outof_15000.rds') 
converged_genes_3000 = names(is_model_converged)[is_model_converged]
Nonconverged_genes_3000 = names(is_model_converged)[!is_model_converged]

conv_df_3000 = rbind(get_convergence_df(converged_genes_3000, UMIs_3000, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_3000, UMIs_3000, tag='non_conv'))

plots_3000 = get_plots(conv_df_3000, run_info = '3000 cells - strain as random effect',
                       is_model_converged =is_model_converged, 
                       perc_non_conv=perc_non_conv_3000)
plots_3000[[6]]
rm(is_model_converged)
######### sample_3000_sampRE_strainRE

is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainRE/Rat_is_model_converged_samleSize_3000_UMIs_sample_RE_strain_RE.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_3000_sampleRE_strainRE = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

converged_genes_3000_sampleRE_strainRE = names(is_model_converged)[is_model_converged]
Nonconverged_genes_3000_sampleRE_strainRE = names(is_model_converged)[!is_model_converged]

conv_df_3000_sampleRE_strainRE = rbind(get_convergence_df(converged_genes_3000_sampleRE_strainRE, UMIs_3000, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_3000_sampleRE_strainRE, UMIs_3000, tag='non_conv'))


plots_3000_sampleRE_strainRE = get_plots(conv_df_3000_sampleRE_strainRE, 
                                         run_info = '3000 cells - sample and strain both as random effect',
                                         is_model_converged =is_model_converged, 
                                         perc_non_conv=perc_non_conv_3000_sampleRE_strainRE)
plots_3000_sampleRE_strainRE[[6]]


######### sample_3000_sampRE_strain as fixed effect

is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainFE/Rat_is_model_converged_samleSize_3000_UMIs_sample_RE_strain_FE.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_3000_sampleRE_strainFE = round(sum(!is_model_converged)*100/length(is_model_converged),3)

converged_genes_3000_sampleRE_strainFE = names(is_model_converged)[is_model_converged]
Nonconverged_genes_3000_sampleRE_strainFE = names(is_model_converged)[!is_model_converged]

conv_df_3000_sampleRE_strainFE = rbind(get_convergence_df(converged_genes_3000_sampleRE_strainFE, UMIs_3000, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_3000_sampleRE_strainFE, UMIs_3000, tag='non_conv'))

plots_3000_sampleRE_strainFE = get_plots(conv_df_3000_sampleRE_strainFE, 
                                         run_info = '3000 cells - sample as random, strain as fixed effect',
                                         is_model_converged =is_model_converged, 
                                         perc_non_conv=perc_non_conv_3000_sampleRE_strainFE)
plots_3000_sampleRE_strainFE[[6]]


############# 3000 strain as random effect - with self implemented Imputation

is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/Imputation/After_imp_3000_selfImp/Rat_is_model_converged_samleSize_3000_UMIs_strain_RE_AfterImpute_selfImp.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_3000_imp = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

UMIs_3000 = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/selected_umi_3000_outof_15000.rds') 
converged_genes_3000_imp = names(is_model_converged)[is_model_converged]
Nonconverged_genes_3000_imp = names(is_model_converged)[!is_model_converged]

conv_df_3000_imp = rbind(get_convergence_df(converged_genes_3000_imp, UMIs_3000, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_3000_imp, UMIs_3000, tag='non_conv'))

plots_3000_imp = get_plots(conv_df_3000_imp, run_info = '3000 cells - strain as random effect - imputated data',
                       is_model_converged =is_model_converged, 
                       perc_non_conv=perc_non_conv_3000_imp)
plots_3000_imp[[6]]
rm(is_model_converged)
#######################################
### how many models did not converge
ggvenn(list(nonConv_3000=Nonconverged_genes_3000,
            nonConv_15000=Nonconverged_genes_15000))

ggvenn(list(Conv_3000=converged_genes_3000,
            Conv_15000=converged_genes_15000))

ggvenn(list(nonConv_3000=Nonconverged_genes_3000,
            nonConv_3000_sampleRE_strainFE=Nonconverged_genes_3000_sampleRE_strainFE,
            nonConv_3000_sampleRE_strainRE=Nonconverged_genes_3000_sampleRE_strainRE), set_name_size=4)

ggvenn(list(Conv_3000=converged_genes_3000,
            Conv_3000_sampleRE_strainFE=converged_genes_3000_sampleRE_strainFE,
            Conv_3000_sampleRE_strainRE=converged_genes_3000_sampleRE_strainRE), set_name_size=4)

ggvenn(list(nonConv_3000=Nonconverged_genes_3000,
            nonConv_3000_sampleRE_strainFE=Nonconverged_genes_3000_sampleRE_strainFE,
            nonConv_3000_sampleRE_strainRE=Nonconverged_genes_3000_sampleRE_strainRE,
            nonConv_15000=Nonconverged_genes_15000), set_name_size=4)

ggvenn(list(Conv_3000=converged_genes_3000,
            Conv_3000_sampleRE_strainFE=converged_genes_3000_sampleRE_strainFE,
            Conv_3000_sampleRE_strainRE=converged_genes_3000_sampleRE_strainRE,
            Conv_15000=converged_genes_15000), set_name_size=4)



length(Nonconverged_genes_3000)
sum(!Nonconverged_genes_3000 %in% Nonconverged_genes_15000)
sum(Nonconverged_genes_3000 %in% Nonconverged_genes_15000)

sum(!Nonconverged_genes_15000 %in% Nonconverged_genes_3000)
length(Nonconverged_genes_15000)


######### Evaluating the effect of Imputation on convergence #########
###### Importing the data and converting it to seurat object ######
data_before_impute = readRDS( '~/scLMM/ratLiver_LMM_results/Before_impute_subset_RatData.rds') 
data_after_impute = readRDS( '~/scLMM/ratLiver_LMM_results/Imputed_subset_RatData.rds') 

data_before_impute[is.na(data_before_impute)] <- 0
data_before_impute = CreateAssayObject(as.matrix(t(data_after_impute)))
dim(data_before_impute)

data_after_impute[is.na(data_after_impute)] <- 0
data_after_impute = CreateAssayObject(as.matrix(t(data_after_impute)))
dim(data_after_impute)

######################################################
################# Before Imputation
is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/Imputation/Before_imp/Rat_is_model_converged_samleSize_3000_UMIs_strain_RE_BeforeImpute.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_beforeImp = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

converged_genes_beforeImp = names(is_model_converged)[is_model_converged]
Nonconverged_genes_beforeImp = names(is_model_converged)[!is_model_converged]

dim(data_before_impute)
UMIs_beforeImp = colnames(data_before_impute)
conv_df_beforeImp = rbind(get_convergence_df(converged_genes_beforeImp, UMIs_beforeImp, tag='conv'), 
                     get_convergence_df(Nonconverged_genes_beforeImp, UMIs_beforeImp, tag='non_conv'))

plots_beforeImp = get_plots(conv_df_beforeImp, run_info = 'Before Imputation - strain as random effect',
                       is_model_converged =is_model_converged, 
                       perc_non_conv=perc_non_conv_beforeImp)

plots_beforeImp[[6]]


################# After Imputation
is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/Imputation/After_imp/Rat_is_model_converged_samleSize_3000_UMIs_strain_RE_AfterImpute.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
perc_non_conv_afterImp = round(sum(!is_model_converged)*100/length(is_model_converged),3) # 29.50042

converged_genes_afterImp = names(is_model_converged)[is_model_converged]
Nonconverged_genes_afterImp = names(is_model_converged)[!is_model_converged]

dim(data_after_impute)
UMIs_afterImp = colnames(data_after_impute)
conv_df_afterImp = rbind(get_convergence_df(converged_genes_afterImp, UMIs_afterImp, tag='conv'), 
                          get_convergence_df(Nonconverged_genes_afterImp, UMIs_afterImp, tag='non_conv'))

plots_afterImp = get_plots(conv_df_afterImp, run_info = 'After Imputation - strain as random effect',
                            is_model_converged =is_model_converged, 
                            perc_non_conv=perc_non_conv_afterImp)

plots_afterImp[[6]]


ggvenn(list(nonConv_afterImp=Nonconverged_genes_afterImp,
           nonConv_beforeImp=Nonconverged_genes_beforeImp))

ggvenn(list(Conv_afterImp=converged_genes_afterImp,
            Conv_beforeImp=converged_genes_beforeImp))




