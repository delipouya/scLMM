### loading the required libraries
source('Codes/Functions.R')
Initialize()

#### new syncronized thresholds: 
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250
# if sample specific QC: DA1: 6623, DA2: 7112 , LEW1: 5457, LEW2: 3844
# if the same QC for all of them: DA1: 6645, DA2: 39799, LEW1: 12124,  LEW2: 7857

# ## import the data
seur_DA_10WK_03 <- CreateSeuratObject(counts=Read10X('Data/rat_DA_M_10WK_003', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")
seur_DA_01 <- CreateSeuratObject(counts=Read10X('Data/rat_DA_01_reseq', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")
seur_LEW_01 <- CreateSeuratObject(counts=Read10X('Data/rat_Lew_01', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")
seur_LEW_02 <- CreateSeuratObject(counts=Read10X('Data/rat_Lew_02', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")

sample_list = list( 'DA_01'= seur_DA_01, 
                    'DA_10WK_03'=seur_DA_10WK_03,
                    'LEW_01'= seur_LEW_01, 
                    'LEW_02' = seur_LEW_02)


MIT_PATTERN = '^Mt-'


sample_list_filt = sapply(1:length(sample_list), function(i){
  
  a_sample = sample_list[[i]]
  INPUT_NAME = names(sample_list)[i]
  
  mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(a_sample))
  a_sample[["mito_perc"]] <- PercentageFeatureSet(a_sample, features = mito_genes_index)
  
  to_drop_mito <- a_sample$mito_perc > MIT_CUT_OFF
  to_drop_lib_size <- a_sample$nCount_RNA < LIB_SIZE_CUT_OFF 
  #to_drop_det_genes <- a_sample$nFeature_RNA < NUM_GENES_DETECTED 
  
  a_sample <- a_sample[,!to_drop_mito & !to_drop_lib_size]
  
  
  return(a_sample)
  
}, simplify = F)

(lapply(sample_list_filt, ncol))


sample_list_filt = sapply(1:length(sample_list), function(i){
  
  a_sample = sample_list[[i]]
  INPUT_NAME = names(sample_list)[i]
  
  if(INPUT_NAME == 'DA_01') {LIB_SIZE_CUT_OFF=1500; MIT_CUT_OFF=30}
  if(INPUT_NAME == 'DA_10WK_03') {LIB_SIZE_CUT_OFF=2000; MIT_CUT_OFF=20}
  if(INPUT_NAME %in% c('LEW_01', 'LEW_02')) {LIB_SIZE_CUT_OFF=2000; MIT_CUT_OFF=40}
  
  mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(a_sample))
  a_sample[["mito_perc"]] <- PercentageFeatureSet(a_sample, features = mito_genes_index)
  
  to_drop_mito <- a_sample$mito_perc > MIT_CUT_OFF
  to_drop_lib_size <- a_sample$nCount_RNA < LIB_SIZE_CUT_OFF 
  to_drop_det_genes <- a_sample$nFeature_RNA < NUM_GENES_DETECTED 
  
  a_sample <- a_sample[,!to_drop_mito & !to_drop_lib_size]
  
  
  return(a_sample)
  
}, simplify = F)


names(sample_list_filt) = names(sample_list)
lapply(sample_list_filt, dim)

merged_samples_count <- merge(sample_list_filt[[1]], c(sample_list_filt[[2]], sample_list_filt[[3]], sample_list_filt[[4]] ), # 
                            add.cell.ids = names(sample_list_filt), 
                            project = "rat_data", 
                            merge.data = TRUE)

colnames(merged_samples_count)
merged_samples_count = readRDS('~/scLMM/set1_merged_samples_count.rds')

dat <- readRDS(file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.rds")
dat <- CreateSeuratObject(dat)
sample_info = sapply(strsplit(colnames(dat), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
strain_info = sapply(strsplit(colnames(dat), '_'), '[[', 1)
table(sample_info)
dat$sample = sample_info
dat$strain = strain_info


library(SeuratDisk)
SaveH5Seurat(dat, filename = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.h5seurat", overwrite = TRUE)
source_file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.h5seurat"
dest_file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData.h5ad"
Convert(source_file, dest_file, assay="RNA", overwrite = TRUE)

