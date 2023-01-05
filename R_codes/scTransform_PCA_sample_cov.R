library(Seurat)
library(ggplot2)
library(sctransform)

rat_h5seur = LoadH5Seurat("~/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5seurat")
pbmc <- rat_h5seur
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = c("sample", 'strain'), verbose = TRUE)
pbmc <- RunPCA(pbmc, verbose = TRUE)
pca_df = data.frame(Embeddings(pbmc, 'pca'))
pca_df$strain = pbmc$strain
pca_df$sample = pbmc$sample
pca_df$cluster = pbmc$cluster
pca_df[1:5,]

ggplot(pca_df, aes(x=PC_1, y=PC_2, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_2, color=sample))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_2, color=strain))+geom_point(size=1)+theme_classic()

ggplot(pca_df, aes(x=PC_1, y=PC_3, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_4, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_5, color=cluster))+geom_point(size=1)+theme_classic()

ggplot(pca_df, aes(x=PC_1, y=PC_6, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_7, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_8, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_9, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_10, color=cluster))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_11, color=strain))+geom_point(size=1)+theme_classic()
ggplot(pca_df, aes(x=PC_1, y=PC_14, color=cluster))+geom_point(size=1)+theme_classic()




#BiocManager::install("glmGamPoi")
#pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc[['pca']]


pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()
