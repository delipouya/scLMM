source('~/RatLiver/Codes/Functions.R')
#source('~/RatLiverCodes/FactorAnalysisUtils.R')
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


Fit_GLMM <- function(i, design_mat, Y){
  
  y = Y[, genelist[i]]
  glmer_data = data.frame(y=y, rand=design_mat$strainLEW)
  formula = 'y ~ 1 + (1 | rand)'
  
  res_glmer <- glmer(formula = as.formula(formula), 
                     data = glmer_data, 
                     family = family)
  return(res_glmer)
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

#merged_samples_count = readRDS('~/scLMM/set1_merged_samples_count.rds')

design_mat = readRDS('~/scLMM/input_data_designMat/designMatrix_rat_set1_countData_refined.rds')
design_mat <- data.frame(design_mat[,1:2])
colnames(design_mat) = c('Intercept', 'strainLEW')

x <- as.matrix(design_mat[,1])	
z <- as.matrix(design_mat$strainLEW) #strain as a random effect
head(z)
table(z)

############## raw counts
input_data = readRDS('~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds')
dim(input_data) #[1] 32883 23036
input_data[1:6, 1:3]
ncells <- rowSums(input_data > 0) ## number of cells which express the gene of interest
hist(ncells)
hist(log2(ncells))
sum(ncells >= 2^4) #[1] 12437
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))

########### filtering gene to include the ones which have been expressed in at least 20 cells
y <- t(input_data[rowSums(input_data > 0) >= 20, ]) 
Y <- as.matrix(y)
genelist = colnames(Y)
family = "poisson"

genelist_org = genelist
genelist = genelist[1:10]
##################
results_list = list()
t1 <- Sys.time()
for (i in 1:length(genelist)){
  
  y = Y[, genelist[i]]
  glmer_data = data.frame(y=y, rand=design_mat$strainLEW)
  formula = 'y ~ 1 + (1 | rand)'
  head(glmer_data)
  
  res_glmer <- glmer(formula = as.formula(formula), 
                     data = glmer_data, 
                     family = family)
  
  summary(res_glmer)
  results_list[[genelist[i]]] = res_glmer
}

t2 <- Sys.time()





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

t2 <- Sys.time()
print(paste0('total time elapsed: ', t2-t1 ))
saveRDS(results_list, '~/scLMM/ratLiver_LMM_results/GLMM_results_lme4_countsData.rds')


#### results interpretation
#v <- vcov.merMod(fitglmer, corr=TRUE)
#as(v, "corMatrix")
#fixef(fitglmer)
#as.data.frame(VarCorr(fitglmer))[, "vcov"]

#####################
### output interpretation
#summary(fitglmer)
# v <- vcov.merMod(res_glmer, corr=TRUE)
#as(v, "corMatrix")
# fixef(res_glmer)
#as.data.frame(VarCorr(fitglmer))[, "vcov"]




### defining the formula in an automatic way 
md <- as.formula(paste0("y ~ 0 + ", 
                        paste0(colnames(X), collapse = "+"), " + (0 + ",
                        paste0(colnames(G), collapse = "+"), " || group)"))

### if no fixed effect included:
md <- as.formula(paste0("y ~ 0 + (0 + ", paste0(colnames(G), collapse = "+"), " || group)"))



