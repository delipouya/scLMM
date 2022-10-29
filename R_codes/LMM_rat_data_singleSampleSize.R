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


Fit_lmer_model <- function(i, merged_samples, data){
  input.data = data.frame(gene=data[i,], 
                          cluster=as.factor(merged_samples$cluster), #[selected_umi_index]
                          sample=as.factor(merged_samples$sample_name),
                          strain=as.factor(merged_samples$strain),
                          numFeatures=merged_samples$nFeature_RNA,
                          librarySize=merged_samples$nCount_RNA
  )
  Model = lmer(gene ~ (1|strain) , data=input.data, REML=T) #sample +(1|sample) + 
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

sample_size = 3000
test_info = ''
# add_info = '' 
# add_info = '_sample_RE_strain_RE'
# add_info = '_sample_RE_strain_FE'
add_info = '_strain_RE_AfterImpute_selfImp'

#### making sure the 3000 subsample is a part of the 15000 sample set 
# umi_samples_15000 = readRDS('ratLiver_LMM_results/subsample_results/sample_15000/Rat_subsampleUMIs_samleSize_150002.rds')
umi_samples_3000 = readRDS('~/scLMM/ratLiver_LMM_results/subsample_results/selected_umi_3000_outof_15000.rds')
merged_samples_sub = merged_samples[,umi_samples_3000]

##########  impute the single cell data - manual implementation ################
data = GetAssayData(merged_samples_sub)
min_data = min(as.matrix(data[data!=0]))
max_data = max(as.matrix(data[data!=0]))
data[data==0] = runif(n=sum(data==0), min = min_data, max = max_data)
sum(data==0)

################################################################
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

##### Open the log #####
log_open(paste0("~/scLMM/ratLiver_LMM_results/subsample_results/rat_subsample_LMMs_",
                sample_size,"_UMIs",add_info, test_info,".log"))

log_print(paste0('sample size: ', sample_size), hide_notes=T)

#### inside loop to fit a LMM model for each gene
log_print(paste0('initial used memory: ', mem_used()), hide_notes=T) 
log_print('Fitting one LMM for each gene started...', hide_notes=T)

iter_to_print = 10
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
  
  
  #if(i %% iter_to_print ==0) {
  #  log_print(paste0('In loop used memory: ', mem_used() ,' - #iteration (gene/trained model index): ', i ), hide_notes=T, blank_after=F)} 
  
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

#saveRDS(models,paste0('~/scLMM/ratLiver_LMM_results/subsample_results/Rat_modelsList_samleSize_', sample_size, "_UMIs",add_info, '.rds')) 
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
log_close()



res <- residuals(mod.lmer)
fit <- fitted(mod.lmer)
