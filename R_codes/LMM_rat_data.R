is_test = F

source('~/RatLiver/Codes/Functions.R')
source('~/RatLiverCodes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(lme4) # load library
library(arm) # convenience functions for regression in R
library(dplyr)
library("pryr")
library(logr)

Fit_lmer_model <- function(i, merged_samples, data){
  input.data = data.frame(gene=data[i,], 
                          cluster=as.factor(merged_samples$cluster), #[selected_umi_index]
                          sample=as.factor(merged_samples$sample_name),
                          strain=as.factor(merged_samples$strain),
                          numFeatures=merged_samples$nFeature_RNA,
                          librarySize=merged_samples$nCount_RNA
  )
  Model = lmer(gene ~ (1|strain) + (1|sample), data=input.data, REML=T) 
  return(Model)
}
get_RandEff_Var <- function(model){
  sum = summary(model)
  RandomEffVar = data.frame(sum$varcor) # sdcor: SD - vcov: Variance
  return(RandomEffVar$vcov[1]) 
}

get_Residual_Var <- function(model){
  sum = summary(model)
  RandomEffVar = data.frame(sum$varcor) # sdcor: SD - vcov: Variance
  return(RandomEffVar$vcov[2]) 
}

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
# feature_size = 500
# feature_size = nrow(merged_samples)
# selected_gene_index = sample(size = feature_size, 1:nrow(merged_samples))
# 
# #Itgal_index = which(rownames(merged_samples) == 'Itgal')
# #data_sub = GetAssayData(merged_samples[c(selected_gene_index, Itgal_index),selected_umi_index])
# merged_samples_org = merged_samples
# merged_samples = merged_samples[,merged_samples$cluster %in% c(5, 10)]

iter_to_print = 500

if(is_test){
  iter_to_print = 2
  merged_samples = merged_samples[1:100,]
}


sample_size_vec = c(seq(3000, ncol(merged_samples), by = 3000), ncol(merged_samples))
sample_size_vec = c(3000) ### the memory limit for the number of cells to handle
sample_size = sample_size_vec[1]
  
#### making sure the 3000 subsample is a part of the 15000 sample set 
umi_samples_15000 = readRDS('ratLiver_LMM_results/subsample_results/sample_15000/Rat_subsampleUMIs_samleSize_150002.rds')

##### Open the log #####

#log_open("~/scLMM/ratLiver_LMM_results/subsample_results/rat_subsample_LMMs_15000_UMIs_2.log")
#log_open("~/scLMM/ratLiver_LMM_results/subsample_results/rat_subsample_LMMs_9000_UMIs_test.log")
add_info = '' 
#add_info = '_sample_RE_strain_RE'
test_info = ''
if(is_test) test_info = '_test'
log_open(paste0("~/scLMM/ratLiver_LMM_results/subsample_results/rat_subsample_LMMs_",sample_size,"_UMIs",add_info, test_info,".log"))


for(sample_size in sample_size_vec){
  
  log_print(paste0('sample size: ', sample_size), hide_notes=T)
  
  if(sample_size <  ncol(merged_samples)){
    selected_umi_index = sample(size = sample_size, 1:ncol(merged_samples))
    #selected_umi_index = readRDS(paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_subsampleUMIs_samleSize_', sample_size, '.rds'))
    merged_samples_sub = merged_samples[,selected_umi_index]
    #saveRDS(selected_umi_index, paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_subsampleUMIs_samleSize_', sample_size, '2.rds'))
    
  }else{
    merged_samples_sub = merged_samples
  }
  
  data = GetAssayData(merged_samples_sub)
  
  ### pre-defining the models list
  models = vector(mode = "list", length = nrow(data))
  is_model_converged = vector(mode = "list", length = nrow(data))
  RandEff_Var_list = vector(mode = "list", length = nrow(data))
  Residual_Var_list = vector(mode = "list", length = nrow(data))
  coef_list = vector(mode = "list", length = nrow(data))
  
  names(models) = rownames(data)
  names(is_model_converged) = rownames(data)
  names(RandEff_Var_list) = rownames(data)
  names(Residual_Var_list) = rownames(data)
  names(coef_list) = rownames(data)
  
  #### inside loop to fit a LMM model for each gene
  log_print(paste0('initial used memory: ', mem_used()), hide_notes=T) 
  log_print('Fitting one LMM for each gene started...', hide_notes=T)
  
  start <- Sys.time()
  for(i in 1:nrow(data)){
    #print(i)
    models[[i]] = Fit_lmer_model(i, merged_samples_sub, data)
    
    is_model_converged[[i]] = models[[i]]@optinfo$conv$lme4
    
    if(isEmpty(is_model_converged[i])) {
      
      RandEff_Var_list[[i]] = get_RandEff_Var(models[[i]])
      Residual_Var_list[[i]] = get_Residual_Var(models[[i]])
      coef_list[[i]] = coef(models[[i]])
      }
    
    
    if(i %% iter_to_print ==0) {
      log_print(paste0('In loop used memory: ', mem_used() ,' - #iteration (gene/trained model index): ', i ), hide_notes=T, blank_after=F)} 
    
    memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))
    if(memfree<500){
      log_print(paste0('The amount of available memory is: '), memfree)
      log_print('Break the loop for lack of available memory')
      break
    }
  }
  
  end <- Sys.time()
  total_time <- as.numeric(end - start, units = "mins") 
  
  log_print('Fitting models is finished.', hide_notes=T)
  log_print(paste0('Final used memory: ', mem_used()), hide_notes=T) 
  
  log_print(paste0('total time requied for sample-size = ', sample_size, ' and num-genes = ', nrow(data), ' : '), hide_notes=T)
  log_print(paste0(round(total_time,5), ' minutues (total dataset) = ', round(total_time*60/nrow(data), 5) ,' seconds per gene'), hide_notes=T)
  
  #saveRDS(models,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_modelsList_samleSize_', sample_size, '.rds')) 
  saveRDS(is_model_converged,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_is_model_converged_samleSize_', sample_size,"_UMIs",add_info, test_info, '.rds')) 
  saveRDS(RandEff_Var_list,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_RandEff_Var_list_samleSize_', sample_size,"_UMIs",add_info, test_info, '.rds')) 
  saveRDS(Residual_Var_list,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_Residual_Var_list_samleSize_', sample_size,"_UMIs",add_info, test_info, '.rds')) 
  saveRDS(coef_list,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_coef_list_samleSize_', sample_size,"_UMIs",add_info, test_info, '.rds')) 
  
  
  ###### releasing memory and checking system's available memory ###### 
  #rm(models)
  #rm(merged_samples_sub)
  #rm(data)
  gc()
  
  ################## 
  log_print('########################################', hide_notes=T)
}

log_close()




########### evaluating each model using the saved models object ###########
models = readRDS(paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_modelsList_samleSize_', sample_size, '.rds')) ## check run_3 screen

is_model_converged <- unlist(lapply(models, function(x) isEmpty(unlist((x@optinfo$conv$lme4)))))
is_model_NOT_converged <- unlist(lapply(models, function(x) !isEmpty(unlist((x@optinfo$conv$lme4)))))
genes_not_converged = names(is_model_NOT_converged[is_model_NOT_converged])

RandEff_Var_list = lapply(models, get_RandEff_Var)
Residual_Var_list = lapply(models, get_Residual_Var)

######### evaluating each model using the saved info extracted from the model objects ########
is_model_converged = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainRE/Rat_is_model_converged_samleSize_3000_UMIs_sample_RE_strain_RE.rds') 
is_model_converged = sapply(1:length(is_model_converged), function(i) isEmpty(is_model_converged[i]))
RandEff_Var_list = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainRE/Rat_RandEff_Var_list_samleSize_3000_UMIs_sample_RE_strain_RE.rds') 
Residual_Var_list = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainRE/Rat_Residual_Var_list_samleSize_3000_UMIs_sample_RE_strain_RE.rds') 
coef_list = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/sample_3000_sampRE_strainRE/Rat_coef_list_samleSize_3000_UMIs_sample_RE_strain_RE.rds')

#######################################
### how many models did not converge
sum(!is_model_converged)*100/length(is_model_converged)

### coef_list contains many NULL objects
sum(sapply(1:length(coef_list), function(i)isEmpty(coef_list[[i]][[1]][1])))/length(coef_list)

RandVar_df = data.frame(gene=names(unlist(RandEff_Var_list)), Rand_Var=unlist(RandEff_Var_list))
RandVar_df_ord = RandVar_df[order(RandVar_df$Rand_Var, decreasing = T),]
rownames(RandVar_df_ord) = NULL
gridExtra::grid.table(head(RandVar_df_ord,20))
dev.off()

Residual_Var_df = data.frame(gene=names(unlist(Residual_Var_list)), Resid_Var=unlist(Residual_Var_list))
Residual_Var_df_ord = Residual_Var_df[order(Residual_Var_df$Resid_Var, decreasing = T),]
rownames(Residual_Var_df_ord) = NULL
gridExtra::grid.table(head(Residual_Var_df_ord,20))



lapply(coef_list[1:10], function(x) isEmpty(x[[1]]))

Variance_df = cbind(RandVar_df, Residual_Var_df)[,-1] #, is_model_converged
head(Variance_df)
ggplot(Variance_df, aes(Rand_Var, Resid_Var, color=Rand_Var))+geom_point()+
  theme_classic()+scale_color_viridis(option = 'viridis',direction = -1)
ggplot(Variance_df, aes(Rand_Var, Resid_Var, color=is_model_converged))+
  geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))

Variance_df = cbind(RandVar_df, Residual_Var_df)[,-1]
merged_df = merge(Variance_df, PCA_loadings, by.x='gene', by.y='gene')
dim(merged_df)
merged_df$PC1_loading_abs <- abs(merged_df$PC1_loading)
pheatmap(cor(merged_df[,-1]))

############# Rand_Var and PC loading
ggplot(merged_df, aes(Rand_Var, PC1_loading_abs, color=Resid_Var))+geom_point()+theme_classic()+scale_color_viridis(option = 'magma',direction = 1)
ggplot(merged_df, aes(Rand_Var, PC1_loading, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
ggplot(merged_df, aes(Rand_Var, PC1_loading_abs, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
############# Resid_Var and PC loading
ggplot(merged_df, aes(Resid_Var, PC1_loading, color=Rand_Var))+geom_point()+theme_classic()+scale_color_viridis(option = 'inferno',direction = 1)
ggplot(merged_df, aes(Resid_Var, PC1_loading, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))
ggplot(merged_df, aes(Resid_Var, PC1_loading_abs, color=is_model_converged))+geom_point(alpha=0.5)+theme_classic()+scale_color_manual(values = c('red', 'grey'))

#############
ggplot(merged_df, aes(Rand_Var, Resid_Var , color=PC1_loading_abs))+geom_point()+theme_classic()+scale_color_viridis(option = 'magma',direction = 1)


###########################################################
############## Parsing out the log file ################## 
########################################################
logFile = readLines(con = '~/scLMM/ratLiver_LMM_results/log/rat_subsample_LMMs_test.log')
logFile = readLines(con = '~/scLMM/ratLiver_LMM_results/subsample_results/log/rat_subsample_LMMs_test.log')
logFile = logFile[-c(1:10)]

sample_size_values = as.numeric(gsub('sample size: ', '', logFile[grepl('sample size:', logFile)]))

time_to_fit = logFile[(grep('total time requied ', logFile))+2]
fit_min_total = unlist(lapply(str_split(time_to_fit, ' = '), function(x) str_split(x[1], ' ')[[1]][1]))
fit_sec_per_gene = unlist(lapply(str_split(time_to_fit, ' = '), function(x) str_split(x[2], ' ')[[1]][1]))

time_df = data.frame(fit_min_total=as.double(fit_min_total), 
                     fit_sec_per_gene=as.double(fit_sec_per_gene), 
                     sample_size_values=sample_size_values[-length(sample_size_values)])

ggplot(time_df, aes(x = sample_size_values, y=fit_min_total))+geom_point(size=3)+ 
  theme_classic()+geom_smooth(method = "lm", se = TRUE)+
  xlab('Sample Size (#UMIs)')+ylab('Total amount of time (min) to fit all genes')+
  ggtitle(paste0('Time to fit LMM models on all genes  \nsample: rat liver (#genes: ', 
                 nrow(your_scRNAseq_data_object), 
                 ' #UMIs=', ncol(your_scRNAseq_data_object),
                 ') \nout of RAM on #UMI:18000 #genes:11500'))

ggplot(time_df, aes(x = sample_size_values, y=fit_sec_per_gene))+geom_point(size=3)+
  theme_classic()+geom_smooth(method = "lm", se = TRUE)+
  xlab('Sample Size (#UMIs)')+ylab('amount of time (seconds) to fit a model per genes')+
  ggtitle(paste0('Time to fit LMM models on a single gene  \nsample: rat liver (#genes: ', 
                 nrow(your_scRNAseq_data_object), 
                 ' #UMIs=', ncol(your_scRNAseq_data_object),
                 ') \nout of RAM on #UMI:18000 #genes:11500'))


memory_usage = as.numeric(unlist(lapply(str_split(logFile[grep('In loop used memory', logFile)], '-'), function(x) str_split(x, ':')[[1]][2])))
iteration_numbers_all = as.numeric(unlist(lapply(str_split(logFile[grep('In loop used memory', logFile)], '-'), function(x) str_split(x, ':')[[2]][2])))
sample_size_annotation = Reduce(c, lapply(sample_size_values, rep, 26, simplify = F))
sample_size_annotation = sample_size_annotation[1:length(memory_usage)] ### No info on the last sample sizes (process killed)

mem_df = data.frame(memory_usage, num_fit_genes_iteration=iteration_numbers_all, sample_size_annotation)
ggplot(mem_df, aes(x=sample_size_annotation, y=memory_usage, color=num_fit_genes_iteration))+
  geom_point()+labs(color="Number of genes \n(= # LMM models)")+
  theme_classic()+geom_line(aes(group = num_fit_genes_iteration))+
  xlab('Sample Size (#UMIs)')+ylab('Memory usage (megabytes)')+
  scale_color_viridis(option = "D")+
  ggtitle(paste0('Memory usage to fit LMM models on a increasing number of genes  \nsample: rat liver (#genes: ', 
                 nrow(your_scRNAseq_data_object), 
                 ' #UMIs=', ncol(your_scRNAseq_data_object),
                 ') \nout of RAM on #UMI:18000 #genes:11500'))



########################################################################
######## evaluating if re-running a model would be effective on the convergence ########
########################################################################

length(genes_not_converged) ## 32 genes which did not converge
dim(merged_samples_sub) ##  first 100 genes of the 9000 subsample
gene_converged_later = c()
counter = 0
for(x in 1:100){
  for(gene in genes_not_converged){
    print(gene)
    i = which(rownames(merged_samples_sub)==gene)
    model = Fit_lmer_model(i, merged_samples_sub, data)
    if(isEmpty(unlist((model@optinfo$conv$lme4)))) {
      counter = counter+1
      gene_converged_later = c(gene_converged_later, gene)
      }
  }
}






