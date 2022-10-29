library(devtools)
library(nebula)
library(matrixStats)
library(pheatmap)
######################################################
################## Nebula code - CJ ##################
####################################################
##design matrix
design_mat = readRDS('~/scLMM/input_data_designMat/designMatrix_rat_set1_countData_refined.rds')
str(design_mat)
dim(design_mat) 
head(design_mat)

## use clusters as fixed effects and strain as a random effect just as an example
#x <- as.matrix(design_mat[,c(1, 6:ncol(design_mat))])	 ## intercept + clusters as fixed effects
x <- as.matrix(design_mat[,c(1, 2)])	 ## intercept + strain
# x <- as.matrix(design_mat[,c(1, 2, 6:ncol(design_mat))])	 ## intercept + strain
head(x)
colnames(x)

z <- as.matrix(design_mat[,2]) #strain as a random effect
z <- as.matrix(design_mat[,c(3,4,5)]) #sample IDs as a random effect
# z <- as.matrix(design_mat[,6:ncol(design_mat)])
head(z)

library(Seurat)
##raw counts
dat <- readRDS(file = "input_data_designMat/inputdata_rat_set1_countData.rds")
dim(dat) #[1] 32883 23036
dat[1:6, 1:3]
ncells <- rowSums(dat > 0)
hist(ncells)
hist(log2(ncells))
sum(ncells >= 2^4) #[1] 12437
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))
genes_2keep = ncells >= 20

y <- as.matrix(dat[genes_2keep, ])

dim(y)
dim(x)
length(z)

##Filtering out low-expressed genes  with counts per cell<0.5%, 
##which is defined by the ratio between the total count of the gene and the number of cells. 
cells_2keep = colSums(y)/nrow(y) >= 5e-3 ### basically all of them!
y = y[,cells_2keep]
x = x[cells_2keep,]
z = z[cells_2keep,]

group_id = sapply(1:nrow(z), function(i) ifelse(sum(z[i,])>0, colnames(z)[z[i,]==1], 'sampleDA_01'),simplify = T) ## sample as RE
# group_id = sapply(1:nrow(z), function(i) ifelse(sum(z[i,])>0, colnames(z)[z[i,]==1], 'cluster0'),simplify = T) ## cluster as RE
#group_id = ifelse(z==0, 'strainDA', 'strainLEW') ## cluster as RE

table(group_id)
colSums(z)

##nebula
#pmm <- nebula(y, as.factor(z), pred = x, offset = colSums(y), model = 'PMM')
### id: random effect design matrix
### x: fixed effect design matrix
# pmm <- nebula(count = y, id = as.factor(z), pred = x[,1], model = 'PMM') # PMM1
# base PMM1: strain as random effect – fixed effect: 1 intercept
# PMM2 sample IDs as random effects and strain as fixed effect 
# PMM3: sample IDs as random effects and strain and cluster as fixed effect 
# PMM4: cluster as random effects and strain as fixed effect 
# PMM5: strain as random effects and cluster as fixed effect 

#### using the default grouping
pmm <- nebula(count = y, id = as.factor(group_id), pred = x, model = 'PMM') 


#### using a new grouping system
data_g = group_cell(count=y,id=group_id, pred=x)
pmm = nebula(data_g$count,data_g$id,pred=data_g$pred)


str(pmm)
#saveRDS(pmm, file = "rat_set1-nebula-pmm2.rds")
pmm = readRDS(file = "rat_set1-nebula-pmm4.rds")

head(pmm$summary)
# The estimated coefficient, standard error and p-value for each predictor.
thr = 15
pmm_df = as.data.frame(pmm$summary)

pmm_table = pmm_df[order(pmm_df$logFC_cluster1, decreasing = T),][1:thr,] ## DA: 0 (pos)- LEW: 1 (neg)
cols_to_check = colnames(pmm_table) %in% c('logFC_strainLEW', "logFC_cluster1",  'p_strainLEW','multi_coef', 'gene')
gridExtra::grid.table(pmm_table[,cols_to_check])

pmm_df$multi_coef = abs(pmm_df$logFC_cluster1 * pmm_df$logFC_strainLEW)
pmm_table = pmm_df[order(pmm_df$multi_coef, decreasing = T),][1:thr,]
gridExtra::grid.table(pmm_table[,cols_to_check])

dev.off()

# The estimated cell-level and subject-level overdispersions 
overdis_df = data.frame(gene=pmm$summary$gene, overdispersion=pmm$overdispersion$Cell)
gridExtra::grid.table(head(overdis_df[order(overdis_df$overdispersion, decreasing = T),], 25))


# More information about the convergence of the algorithm for each gene. 
# A value of -20 or -30 indicates a potential failure of the convergence.
table(as.character(pmm$convergence))

pmm$random_effect # The subject-level random effects.

#The covariance matrix for the estimated log(FC).
table(pmm$convergence)

pmm_df$convergence = pmm$convergence
gene_stats = data.frame(gene_name=rownames(dat),
                        total_count=rowSums(dat),
                        mean=rowMeans(as.matrix(dat)),
                        variance=rowVars(as.matrix(dat)))
merged_df = merge(gene_stats, pmm_df, by.x = 'gene_name', by.y='gene', all.x = F, all.y=T)


head(merged_df)
dim(merged_df)
dim(gene_stats)
library(ggplot2)
merged_df2 = merged_df[,!colnames(merged_df)%in% c("gene_name","gene_id")]
head(merged_df2)
pheatmap::pheatmap(cor(merged_df2[,-ncol(merged_df2)]))

ggplot(merged_df2, aes(variance, logFC_strainLEW))+geom_point(size=3, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(mean, logFC_strainLEW))+geom_point(size=3, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(convergence, mean, color=logFC_strainLEW))+geom_point(size=3, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(convergence, variance, color=logFC_strainLEW))+geom_point(size=3, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(convergence, variance, color=logFC_strainLEW))+geom_point(size=3, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(p_strainLEW, variance, color=logFC_strainLEW))+geom_point(size=1, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(p_strainLEW, convergence, color=logFC_strainLEW))+geom_point(size=1, alpha=0.5)+theme_classic()
ggplot(merged_df2, aes(p_strainLEW, convergence, color=logFC_strainLEW))+geom_point(size=1, alpha=0.5)+theme_classic()


var15 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC15_loadings.rnk')
colnames(var15) = c('gene', 'var15_load')
var5 = read.table('~/RatLiver/Results/strain_variation_loadings/ranked_files/old_samples_rot_PC5_loadings.rnk')
colnames(var5) = c('gene', 'var5_load')


logfc_thr_pos = 3
logfc_thr_neg = -5.5
merged_var15 = merge(merged_df, var15, by.x = 'gene_name', by.y='gene')
merged_varimax = merge(merged_var15, var5, by.x = 'gene_name', by.y='gene')
head(merged_varimax)
merged_varimax$convergence = as.character(merged_varimax$convergence)
ggplot(merged_varimax, aes(logFC_strainLEW, var15_load, color=p_strainLEW))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula logFC')
ggplot(merged_varimax, aes(logFC_strainLEW, var15_load, color=convergence))+geom_point(size=1.5, alpha=0.8)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula logFC')+
  geom_text(aes(label=ifelse(var15_load>0.1 | var15_load<(-0.1)|logFC_strainLEW>logfc_thr_pos|logFC_strainLEW<(logfc_thr_neg),
                             as.character(gene_name),'')),hjust=0,vjust=0)
ggplot(merged_varimax, aes(p_strainLEW, var15_load, color=convergence))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula p-value')+
  geom_text(aes(label=ifelse(var15_load>0.1|var15_load<(-0.1) ,as.character(gene_name),'')),hjust=0,vjust=0)


ggplot(merged_varimax, aes(logFC_strainLEW, var5_load, color=p_strainLEW))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula logFC')
ggplot(merged_varimax, aes(logFC_strainLEW, var5_load, color=convergence))+geom_point(size=1.5, alpha=0.8)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula logFC')+
  geom_text(aes(label=ifelse(var5_load>0.09 | var5_load<(-0.09)|logFC_strainLEW>logfc_thr_pos|logFC_strainLEW<(logfc_thr_neg),as.character(gene_name),'')),hjust=0,vjust=0)
ggplot(merged_varimax, aes(p_strainLEW, var5_load, color=convergence))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula p-value')+
  geom_text(aes(label=ifelse(var5_load>0.1|var5_load<(-0.1) ,as.character(gene_name),'')),hjust=0,vjust=0)



#### comparing overdispersion rate of model-1 with model3
pmm1 = readRDS(file = "rat_set1-nebula-pmm1.rds")
pmm2 = readRDS(file = "rat_set1-nebula-pmm2.rds")

overdis_df_pmm1 = data.frame(gene=pmm1$summary$gene, overdispersion=pmm1$overdispersion, logfc_1=pmm1$overdispersion)
logFC_pmm2 = data.frame(gene=pmm2$summary$gene, logFC=pmm2$summary$logFC_strainLEW, abs_logFC=abs(pmm2$summary$logFC_strainLEW))
merged_pmm1_2 = merge(overdis_df_pmm1, logFC_pmm2, by.x='gene', by.y='gene')
head(merged_pmm1_2)
ggplot(merged_pmm1_2, aes(x=overdispersion, y=logFC, color=logfc_1))+geom_point()+
  theme_classic()+xlab('gene overdispersion of base model - strain as RE')+
  ylab('strain fixed effect LogFC of model-2 - strain as FE - sample as RE')+
  geom_text(aes(label=ifelse(overdispersion>1 ,as.character(gene),'')),hjust=0,vjust=0, size=3)

ggplot(merged_pmm1_2, aes(x=overdispersion, y=abs_logFC, color=logfc_1))+geom_point()+
  theme_classic()+xlab('gene overdispersion of base model - strain as RE')+
  ylab('strain fixed effect LogFC of model-2 - strain as FE - sample as RE')+
  geom_text(aes(label=ifelse(overdispersion>1 ,as.character(gene),'')),hjust=0,vjust=0, size=3)

## PMM1: base model - Base GLMM: strain as random effect – fixed effect: 1 intercept 
##  nebula(count = y, id = as.factor(z), pred = x[,1], model = 'PMM’)


