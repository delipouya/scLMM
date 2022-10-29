#### Check the tutorial on the package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splat_params.html
library(splatter)
library(scater)
library(purrr) # v. 0.3.4
library(broom) # v. 0.5.6
library(dplyr) # v. 1.0.0
library(ggplot2) # v. 3.3.1
library(Seurat)
library(SingleCellExperiment)
set.seed(1)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')

################################################################
############## Simulating data using the splatter library ###### 
################################################################
#### Simulating data based on leaned parameters from healthy rat liver dataset #####
counts <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
sce <- SingleCellExperiment(list(counts=as.matrix(counts)))
row_sample_num = round(nrow(sce)/10)
col_sample_num = round(ncol(sce)/10)
sce <- sce[sample(1:nrow(sce), row_sample_num, replace=F), 
           sample(1:ncol(sce), col_sample_num, replace=F)]
dim(sce)
# Estimate parameters from rat data
params_rat <- splatEstimate(sce)
params = params_rat

######### Simulating data based on a mock dataset #####
# Create mock data
sce_mock <- mockSCE()
# Estimate parameters from mock data
params_mock <- splatEstimate(sce_mock)
#> NOTE: Library sizes have been found to be normally distributed instead of log-normal. You may want to check this is correct.
# Simulate data using estimated parameters
params <- params_mock
####################################################

#params <- newSplatParams()
params <- setParam(params, "nGenes", 2000)

total_num_cells = 5000
batch1 = 3000
group1_p = 0.3
group2_p = 0.5
batch_factor = 0.08


sim <- splatSimulate(params, batchCells = c(batch1, total_num_cells-batch1), 
                             group.prob = c(group1_p, group2_p, 1-(group1_p+group2_p)),
                     batch.facLoc = batch_factor, batch.facScale = batch_factor, 
                     method = "groups", verbose = T)
dim(sim)
sim <- logNormCounts(sim)
sim <- runPCA(sim)

plotPCA(sim, shape_by = "Batch", colour_by = "Group")
plotPCA(sim, shape_by = "Group", colour_by = "Batch")

sim <- runUMAP(sim)
plotUMAP(sim, shape_by = "Batch", colour_by = "Group")
plotUMAP(sim, shape_by = "Group", colour_by = "Batch")

#BatchFac[Batch] The batch effects factor for each gene for a particular batch.
#DEFac[Group] The differential expression factor for each gene in a particular group. Values of 1 indicate the gene is not differentially expressed.
#SigmaFac[Path] Factor applied to genes that have non-linear changes in expression along a path.

counts(sim)[1:5, 1:5] ## count data is simulated
head(rowData(sim)) ### can be used to compare with coefficients and the residual based factors
head(colData(sim)) 
names(assays(sim))
assays(sim)$CellMeans[1:5, 1:5]

### for models using this data and use batch and group for both residual-based control and covariate 
saveRDS(sim, 'simulated_data/simulated_data_3groups_2batch.rds')




##### simulating from base R functions
ngroup = 2
nrep = 10
b0 = 5
b1 = -2
sd = 2
group = rep( c("group1", "group2"), each = nrep) 
eps = rnorm(n = ngroup*nrep, mean = 0, sd = sd) 
growth = b0 + b1*(group == "group2") + eps 
dat = data.frame(group, growth)
growthfit = lm(growth ~ group, data = dat)
summary(growthfit)


