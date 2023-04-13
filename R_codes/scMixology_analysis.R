source('~/RatLiver/Codes/Functions.R')
Initialize()


###################################################################
load("data/sincell_with_class.RData") ## 3 cell line data
#sce10x_qc contains the read counts after quality control processing from the 10x platform. 
#sce4_qc contains the read counts after quality control processing from the CEL-seq2 platform. 
#scedrop_qc_qc contains the read counts after quality control proessing from the Drop-seq platform.
data_list_3cl = list(sc_10x=sce_sc_10x_qc, 
                     CELseq2=sce_sc_CELseq2_qc, 
                     Dropseq=sce_sc_Dropseq_qc)


#### using the count data and applying seurat normalization to each sample individually
data_list_3cl = lapply(data_list_3cl, function(x) {
  ### extracting the raw counts
  x_count = as.Seurat(x, data = "counts")
  x_count_data <- GetAssayData(object =  x_count[['originalexp']], slot = 'data')
  x_count[["RNA"]] <- CreateAssayObject(data = x_count_data )
  
  ### extracting the log normliazed counts
  x_lognorm = as.Seurat(x, data = "logcounts")
  x_lognorm_data <- GetAssayData(object =  x_lognorm[['originalexp']], slot = 'data')
  #x_count[["lognorm"]] <- CreateAssayObject(data = x_lognorm_data )
  x_count[["SCT"]] <- CreateAssayObject(data = x_lognorm_data )
  
  #### setting the default assay as the RNA slot
  DefaultAssay(x_count) <- "RNA"
  x_count[['originalexp']] <- NULL
  
  ### sctransform does not work - not sure why
  # x = SCTransform(x, conserve.memory=F, verbose=T, return.only.var.genes=F, variable.features.n = nrow(x), 
  #                do.scale = FALSE, do.center = TRUE) ### default value 
  
  return(x_count)
  })


names(data_list_3cl) = c('sc_10X', 'CELseq2', 'Dropseq')
data_list_3cl = sapply(1:length(data_list_3cl), 
                       function(i) {data_list_3cl[[i]]$sample=names(data_list_3cl)[i]; data_list_3cl[[i]]}, simplify = F)
names(data_list_3cl) = c('sc_10X', 'CELseq2', 'Dropseq')


scMix_3cl_merged <- merge(data_list_3cl[[1]], c(data_list_3cl[[2]], data_list_3cl[[3]]),
                          add.cell.ids = names(data_list_3cl), 
                          project = "scMix_3cl", 
                          merge.data = TRUE)

dim(scMix_3cl_merged)
scMix_3cl_merged_raw = GetAssayData(scMix_3cl_merged, assay = 'RNA')
scMix_3cl_merged_indiv.lognorm = GetAssayData(scMix_3cl_merged, assay = 'SCT')


sum(colSums(scMix_3cl_merged_raw) == 0)
sum(colSums(scMix_3cl_merged_indiv.lognorm) == 0)

df = data.frame(colsum_raw=colSums(scMix_3cl_merged_raw),
                colsum_lognorm_indiv=colSums(scMix_3cl_merged_indiv.lognorm),
                colsum_norm_merge=colSums(NormalizeData(scMix_3cl_merged, verbose = TRUE) ))

DefaultAssay(scMix_3cl_merged) <- 'SCT'
df$colsum_scale_merge=colSums(ScaleData(scMix_3cl_merged, verbose = TRUE))
summary(df) #### scaling and not scaling seem to have no effect on the data


variable_features = lapply(data_list_3cl, function(x) VariableFeatures(FindVariableFeatures(x)))
sum(!unique(unlist(variable_features)) %in% rownames(scMix_3cl_merged))
VariableFeatures(scMix_3cl_merged) = unique(unlist(variable_features))


scMix_3cl_merged = ScaleData(scMix_3cl_merged, verbose = TRUE)


scMix_3cl_merged <- RunPCA(scMix_3cl_merged,verbose=T)
plot(100 * scMix_3cl_merged@reductions$pca@stdev^2 / scMix_3cl_merged@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

num_PCs = 40
#### UMAP on PC components and clustering
scMix_3cl_merged <- RunUMAP(scMix_3cl_merged, reduction = "pca", dims = 1:num_PCs, verbose = TRUE)

table(scMix_3cl_merged$sample)
#### UMAP on corrected PC components and clustering
scMix_3cl_merged <- RunHarmony(scMix_3cl_merged, group.by.vars = "sample")
scMix_3cl_merged <- RunUMAP(scMix_3cl_merged, reduction = "harmony", dims = 1:num_PCs, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:num_PCs, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)




df_umap <- data.frame(UMAP_1=getEmb(scMix_3cl_merged, 'umap_h')[,1], 
                      UMAP_2=getEmb(scMix_3cl_merged, 'umap_h')[,2], 
                      library_size= scMix_3cl_merged$nCount_originalexp, 
                      n_expressed=scMix_3cl_merged$nFeature_originalexp ,
                      cell_line=scMix_3cl_merged$cell_line,
                      clusters= as.character(scMix_3cl_merged$RNA_snn_res.0.6), 
                      sample_name=scMix_3cl_merged$sample, 
                      cell_line_demuxlet=scMix_3cl_merged$cell_line_demuxlet
)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(alpha=0.7)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_line_demuxlet))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_line))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)



saveRDS(scMix_3cl_merged, "~/scLMM/sc_mixology/scMix_3cl_merged_analyzed.rds")
###################################################################




load("data/sincell_with_class_5cl.RData") ## 5 cell line data
data_list_5cl = list(sc_10x=sce_sc_10x_5cl_qc, 
                     CELseq2_p1=sc_Celseq2_5cl_p1, 
                     CELseq2_p2=sc_Celseq2_5cl_p2,
                     CELseq2_p3=sc_Celseq2_5cl_p3)

#### using the count data and applying seurat normalization to each sample individually
data_list_5cl = lapply(data_list_5cl, function(x) {
  ### extracting the raw counts
  x_count = as.Seurat(x, data = "counts")
  x_count_data <- GetAssayData(object =  x_count[['originalexp']], slot = 'data')
  x_count[["RNA"]] <- CreateAssayObject(data = x_count_data )
  
  ### extracting the log normliazed counts
  
  x_lognorm = logNormCounts(x)
  x_lognorm_data <- assay(x_lognorm, 'logcounts')
  #x_count[["lognorm"]] <- CreateAssayObject(data = x_lognorm_data )
  x_count[["SCT"]] <- CreateAssayObject(data = x_lognorm_data )
  
  #### setting the default assay as the RNA slot
  DefaultAssay(x_count) <- "RNA"
  x_count[['originalexp']] <- NULL
  
  ### sctransform does not work - not sure why
  # x = SCTransform(x, conserve.memory=F, verbose=T, return.only.var.genes=F, variable.features.n = nrow(x), 
  #                do.scale = FALSE, do.center = TRUE) ### default value 
  
  return(x_count)
})



data_list_5cl_names = names(data_list_5cl)

data_list_5cl = sapply(1:length(data_list_5cl), 
                       function(i) {data_list_5cl[[i]]$sample=names(data_list_5cl)[i]; data_list_5cl[[i]]}, simplify = F)
names(data_list_5cl) = data_list_5cl_names



scMix_5cl_merged <- merge(data_list_5cl[[1]], c(data_list_5cl[[2]], data_list_5cl[[3]], data_list_5cl[[4]]),
                          add.cell.ids = names(data_list_5cl), 
                          project = "scMix_5cl", 
                          merge.data = TRUE)

dim(scMix_5cl_merged)
scMix_5cl_merged_raw = GetAssayData(scMix_5cl_merged, assay = 'RNA')
scMix_5cl_merged_indiv.lognorm = GetAssayData(scMix_5cl_merged, assay = 'SCT')

sum(colSums(scMix_5cl_merged_raw) == 0)
sum(colSums(scMix_5cl_merged_indiv.lognorm) == 0)

df = data.frame(colsum_raw=colSums(scMix_5cl_merged_raw),
                colsum_lognorm_indiv=colSums(scMix_5cl_merged_indiv.lognorm),
                colsum_norm_merge=colSums(NormalizeData(scMix_5cl_merged, verbose = TRUE) ))

DefaultAssay(scMix_5cl_merged) <- 'SCT'
df$colsum_scale_merge=colSums(ScaleData(scMix_5cl_merged, verbose = TRUE))
summary(df) #### scaling and not scaling seem to have no effect on the data


variable_features = lapply(data_list_5cl, function(x) VariableFeatures(FindVariableFeatures(x)))
sum(!unique(unlist(variable_features)) %in% rownames(scMix_5cl_merged))
VariableFeatures(scMix_5cl_merged) = unique(unlist(variable_features))


scMix_5cl_merged = ScaleData(scMix_5cl_merged, verbose = TRUE)


scMix_5cl_merged <- RunPCA(scMix_5cl_merged,verbose=T)
plot(100 * scMix_5cl_merged@reductions$pca@stdev^2 / scMix_5cl_merged@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

num_PCs = 30
#### UMAP on PC components and clustering
scMix_5cl_merged <- RunUMAP(scMix_5cl_merged, reduction = "pca", dims = 1:num_PCs, verbose = TRUE)

table(scMix_5cl_merged$sample)
str_split(scMix_5cl_merged$sample, pattern = '_')

scMix_5cl_merged$protocol = ifelse(scMix_5cl_merged$sample=='sc_10x', 'sc_10x', 'CELseq2')
table(scMix_5cl_merged$protocol, scMix_5cl_merged$sample)

#### UMAP on corrected PC components and clustering
scMix_5cl_merged <- RunHarmony(scMix_5cl_merged, group.by.vars = "sample")
scMix_5cl_merged <- RunUMAP(scMix_5cl_merged, reduction = "harmony", dims = 1:num_PCs, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:num_PCs, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)


df_umap <- data.frame(UMAP_1=getEmb(scMix_5cl_merged, 'umap_h')[,1], 
                      UMAP_2=getEmb(scMix_5cl_merged, 'umap_h')[,2], 
                      library_size= scMix_5cl_merged$nCount_originalexp, 
                      n_expressed=scMix_5cl_merged$nFeature_originalexp ,
                      cell_line=scMix_5cl_merged$cell_line,
                      clusters= as.character(scMix_5cl_merged$SCT_snn_res.0.6), 
                      sample_name=scMix_5cl_merged$sample, 
                      cell_line_demuxlet=scMix_5cl_merged$cell_line_demuxlet
)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(size=1,alpha=0.7)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_line_demuxlet))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_line))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)

saveRDS(scMix_5cl_merged, "~/scLMM/sc_mixology/scMix_5cl_merged_analyzed.rds")


