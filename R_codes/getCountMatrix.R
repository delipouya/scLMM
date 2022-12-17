### loading the required libraries
source('~/RatLiver/Codes/Functions.R')
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


##### adding the cluster information to the metadata
library(Seurat)
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
rm(your_scRNAseq_data_object); gc()
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'LEW-1', 'LEW-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
processed_df = data.frame(cluster=merged_samples$cluster, 
                          sample_name=merged_samples$sample_name, 
                          strain=merged_samples$strain, 
                          cell_id=colnames(merged_samples))
###### refining the cell ID of the processed data ######
sample_info = sapply(strsplit(processed_df$cell_id, '_'), function(x) paste0(x[2], '_' ,x[3]))
pure_cell_id = sapply(strsplit(processed_df$cell_id, '_'), function(x) paste0(x[length(x)]))
sample_info[sample_info=='DA_M'] = 'DA_02'
sample_info[sample_info=='Lew_01'] = 'LEW_01'
sample_info[sample_info=='Lew_02'] = 'LEW_02'
table(sample_info)
processed_df$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
dim(processed_df)

##### refining the cell IDs of the raw data ######
dat_metadata = data.frame(dat@meta.data)
head(dat_metadata)
sum(!rownames(dat_metadata) %in% processed_df$cell_id_2)

sample_info = sapply(strsplit(rownames(dat_metadata), '_'), function(x) paste0(x[1], '_' ,x[2]))
sample_info[sample_info=='DA_10WK'] = 'DA_02'
table(sample_info)
pure_cell_id = sapply(strsplit(rownames(dat_metadata), '_'), function(x) paste0(x[length(x)]))
dat_metadata$cell_id_2 = paste0(sample_info, '_', pure_cell_id)
sum(!dat_metadata$cell_id_2 %in% processed_df$cell_id_2)

merged_df = merge(dat_metadata, processed_df, by.x='cell_id_2', by.y='cell_id_2')
head(merged_df)
sum(merged_df$cell_id_2 != dat_metadata$cell_id_2) ### the order of rows in the merged dataframe is conserved 
### add the cluster information to the raw data's metadata dataframe
dat_metadata$cluster = merged_df$cluster
dat$cluster  = dat_metadata$cluster
dat$refined_cell_ID = dat_metadata$cell_id_2
  
  
library(SeuratDisk)
SaveH5Seurat(dat, filename = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5seurat", overwrite = TRUE)
source_file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5seurat"
dest_file = "~/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5ad"
Convert(source_file, dest_file, assay="RNA", overwrite = TRUE)

