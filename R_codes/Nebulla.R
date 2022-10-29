library(devtools)
library(nebula)
library(matrixStats)
library(pheatmap)
input_data = readRDS('input_data_designMat/inputdata_rat_set1_countData.rds')
designMatrix = readRDS('input_data_designMat/designMatrix_rat_set1_countData.rds')


cells_to_include = colSums(input_data)!=0
genes_to_include = rowSums(input_data)!=0 &  rowVars(as.matrix(input_data))!=0
input_data2 = input_data[genes_to_include, cells_to_include]
designMatrix2 = designMatrix[cells_to_include,]



num_cells_to_include = nrow(designMatrix2)
re = nebula(count = input_data2[,1:num_cells_to_include],
            id = designMatrix2$strainLEW,
            pred=designMatrix2[1:num_cells_to_include,][,c(1,2)])


head(re$summary)
# The estimated coefficient, standard error and p-value for each predictor.
nebula_result_df = re$summary

DA_table = nebula_result_df[order(nebula_result_df$logFC_strainLEW, decreasing = F),][1:15,]
LEW_table = nebula_result_df[order(nebula_result_df$logFC_strainLEW, decreasing = T),][1:15,]

gridExtra::grid.table(LEW_table)
# The estimated cell-level and subject-level overdispersions 
head(re$overdispersion)

# More information about the convergence of the algorithm for each gene. 
# A value of -20 or -30 indicates a potential failure of the convergence.
table(re$convergence)

re$random_effect # The subject-level random effects.

#The covariance matrix for the estimated log(FC).
re$covariance

nebula_result_df$convergence = re$convergence
gene_stats = data.frame(gene_name=rownames(input_data2),
                        total_count=rowSums(input_data2),
                        mean=rowMeans(as.matrix(input_data2)),
                        variance=rowVars(as.matrix(input_data2)))
merged_df = merge(gene_stats, nebula_result_df, by.x = 'gene_name', by.y='gene', all.x = F, all.y=T)

head(merged_df)
dim(merged_df)
dim(gene_stats)
library(ggplot2)
merged_df2 = merged_df[,!colnames(merged_df)%in% c("gene_name","gene_id")]
head(merged_df2)
pheatmap::pheatmap(cor(merged_df2)[1:3,4:10])

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

merged_var15 = merge(merged_df, var15, by.x = 'gene_name', by.y='gene')
merged_varimax = merge(merged_var15, var5, by.x = 'gene_name', by.y='gene')
head(merged_varimax)
merged_varimax$convergence = as.character(merged_varimax$convergence)
ggplot(merged_varimax, aes(logFC_strainLEW, var15_load, color=p_strainLEW))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula logFC')
ggplot(merged_varimax, aes(logFC_strainLEW, var15_load, color=convergence))+geom_point(size=1.5, alpha=0.8)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula logFC')+
  geom_text(aes(label=ifelse(var15_load>0.1 | var15_load<(-0.1)|logFC_strainLEW>10|logFC_strainLEW<(-10) ,as.character(gene_name),'')),hjust=0,vjust=0)
ggplot(merged_varimax, aes(p_strainLEW, var15_load, color=convergence))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 15 correlation with nebula p-value')+
  geom_text(aes(label=ifelse(var15_load>0.1|var15_load<(-0.1) ,as.character(gene_name),'')),hjust=0,vjust=0)


ggplot(merged_varimax, aes(logFC_strainLEW, var5_load, color=p_strainLEW))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula logFC')
ggplot(merged_varimax, aes(logFC_strainLEW, var5_load, color=convergence))+geom_point(size=1.5, alpha=0.8)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula logFC')+
  geom_text(aes(label=ifelse(var5_load>0.09 | var5_load<(-0.09)|logFC_strainLEW>10|logFC_strainLEW<(-10) ,as.character(gene_name),'')),hjust=0,vjust=0)
ggplot(merged_varimax, aes(p_strainLEW, var5_load, color=convergence))+geom_point(size=1.5, alpha=0.5)+
  theme_classic()+ggtitle('varimax 5 correlation with nebula p-value')+
  geom_text(aes(label=ifelse(var5_load>0.1|var5_load<(-0.1) ,as.character(gene_name),'')),hjust=0,vjust=0)



head(input_data)
head(designMatrix)

dim(input_data)
dim(designMatrix)

data(sample_data)
dim(sample_data$count)
sample_data$count[1:5,1:5]
head(sample_data$sid)
length(sample_data$sid)
#> [1] "1" "1" "1" "1" "1" "1"
table(sample_data$sid)
head(sample_data$pred)
dim(sample_data$pred)
df = model.matrix(~X1+X2+cc, data=sample_data$pred)
head(df)



data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)
re = nebula(data_g$count,data_g$id,pred=data_g$pred)

