test = readRDS('~/scLMM/LMM-scRNAseq-jan2024/sclmm_test_kidneyhuman_expcounts.rds')
head(test)

negbn = readRDS('~/scLMM/LMM-scRNAseq-jan2024/nebula_kidneyhuman_expcounts.rds')
negbn <- readRDS('~/scLMM/LMM-scRNAseq-jan2024/nebula_kidneyhuman_expcounts_NBGMM.rds')
##nebula outputs
##summary (statistics): 
##The estimated coefficient, standard error and p-value for each predictor.
str(negbn)

plot(negbn$overdispersion$Subject)
plot(negbn$overdispersion$Cell)
table(negbn$convergence)



##t-values
tvlmm <- test[, grep("_t", colnames(test)), drop = F]
dim(tvlmm)
head(tvlmm)

##p-values
plmm <- test[, grep("_pvalue", colnames(test)), drop = F]
dim(plmm)
head(plmm)


##########
##nebula
rtnebula = NA
rtnebula[[1]] = 10257.0706329346
#Time difference of 3209.861 secs
rtnebula[[1]]/rtlmm[[1]]
#[1] 64.23869
str(negbn)


st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

any(is.na(st))
sum(is.na(st))
st[rowSums(is.na(st)) > 0, ]
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)


##fixed effects, se, and p-values
iFC <- grep("logFC_", colnames(st))
ise <- grep("se_", colnames(st))
ipv <- setdiff(grep("p_", colnames(st)), c(iFC, ise))

##fixed effects
b <- st[, iFC]

##se
se <- st[, ise]

##t-values
tv <- b/se

##p-values
pv <- st[, ipv]

range(pv - 2*pnorm(-abs(tv)), na.rm = T)
#[1] -2.331468e-15  2.386980e-15


##########
##comparison of lmmfit and nebula
##for the genes that lmmfit and nebula are convergent.

all(rownames(tv) == rownames(tvlmm))
any(is.na(tv))
any(is.na(tvlmm))

head(tv)
#### removing Podocytes
tv = tv[,!grepl('Podocyte', colnames(tv))]
tvlmm = tvlmm[,!grepl('Podocyte', colnames(tvlmm))]
pv = pv[,!grepl('Podocyte', colnames(pv))]
plmm = plmm[,!grepl('Podocyte', colnames(plmm))]



hist(negbn$overdispersion$Subject)
hist(negbn$overdispersion$Cell)
##genes with convergence
i <- (negbn$overdispersion$Subject < 4) ## used to be 1 
sum(i)
i <- i & (negbn$overdispersion$Cell < 100) & (negbn$convergence >= -10)
sum(i)
i <- i & (colSums(abs(fit$dlogL) < 1e-5) == nrow(fit$dlogL))
sum(i)
i <- i & (apply(!is.na(tv), 1, all)) & (apply(!is.na(tvlmm), 1, all))
sum(i)



##t-values
j <- 2:ncol(tv)
plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]), 
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1)

##p-values
j <- 2:ncol(pv)
plot(as.matrix(-log10(plmm[i, j])), as.matrix(-log10(pv[i, j])), 
     xlab = "lmmfit -log10(p-values)", ylab = "nebula -log10(p-values)", cex = 0.6)
abline(0, 1)


old.par <- par(no.readonly = TRUE)
on.exit(par(old.par))
##cell-type specific
jset <- grep(":Male", colnames(pv))
nc <- round(sqrt(length(jset)))
nr <- ceiling(length(jset)/nc)
par(mfrow = c(nr, nc), mar = c(5.1,4.1,2.1,1.1))
for (j in jset[1:15]){
  nm <- gsub("p", "+", gsub("_t", "", colnames(tvlmm)[j]))
  plot(tvlmm[i,j], tv[i,j], cex = 0.6,
       xlab = "lmmfit t-values", ylab = "nebula t-values", main = nm, cex.main = 0.8)
  abline(0,1)
}

j = jset[1:15][1]
length(jset)
head(tvlmm)
head(tv)

colnames(tvlmm)[j]
colnames(tv)[j]
nm <- gsub("p", "+", gsub("_t", "", colnames(tvlmm)[j]))
plot(tvlmm[i,j], tv[i,j], cex = 0.6,
     xlab = "lmmfit t-values", ylab = "nebula t-values", main = nm, cex.main = 0.8)

tvlmm[i,j]
tv[i,j]

plot(unlist(tvlmm[i,]), unlist(tv[i,]), cex = 0.6,
     xlab = "lmmfit t-values", ylab = "nebula t-values", main ='all cell types', cex.main = 0.8)

## t-values
j <- 2:ncol(tv) ### exclude intercept
colnames(tv)
plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]), 
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1)

plot(unlist(tvlmm[i, j]), unlist(tv[i, j]), 
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)

library(reshape2)
tvlmm_df = data.frame(melt(tvlmm[i, j]))
colnames(tvlmm_df)[3] = 'sclmm_tval'

tv_df = (as.matrix(tv))
head(tv_df)
dim(tv_df)
tv_df = data.frame(melt(tv_df[i, j]))
head(tv_df)

length(tvlmm_df$value)
length(tv_df$value)
tv_df$Var1 == tvlmm_df$Var1
head(tv_df)
colnames(tv_df)[3] = 'nebula_tval'

merged_df = cbind(tv_df, tvlmm_df)
head(merged_df)
colnames(merged_df)[1] = 'genes'
colnames(merged_df)[2] = 'col_name'
merged_df$nebula_tval = as.numeric(merged_df$nebula_tval)
merged_df$sclmm_tval = as.numeric(merged_df$sclmm_tval)

library(RColorBrewer)
n <- 40
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(merged_df, aes(x=sclmm_tval, y=nebula_tval, color=Var2))+
  geom_point()+theme_classic()+scale_color_manual(values =col_vector )+
  geom_hline(yintercept=5, linetype="dashed", color = "red")+
  geom_vline(xintercept=10, linetype="dashed", color = "red")

df_inconsist = merged_df[merged_df$nebula_tval<1 & merged_df$nebula_tval>(-1) & merged_df$sclmm_tval > 10,] 
head(df_inconsist)
df_inconsist$score = abs(df_inconsist$sclmm_tval/df_inconsist$nebula_tval)
df_inconsist_ord = df_inconsist[order(df_inconsist$score, decreasing = T),]
#write.csv(df_inconsist_ord, '~/scLMM/LMM-scRNAseq-jan2024/inconsistent_nebula_scLMM-kidneyhuman.csv')
write.csv(df_inconsist_ord, '~/scLMM/LMM-scRNAseq-jan2024/inconsistent_nebula_scLMM-kidneyhuman_NBGMM.csv')
ggplot(df_inconsist, aes(x=sclmm_tval, y=nebula_tval, color=Var2))+
  geom_point()+theme_classic()+scale_color_manual(values =col_vector )
ggplot(df_inconsist, aes(x=sclmm_tval, y=nebula_tval, color=score))+
  geom_point()+theme_classic()+scale_color_viridis_b(direction = -1)

df_inconsist_ord

###########################################################################################
############ evaluating the correlation between inconsistency score and variance/sumcounts
##############################################################################
dirData = '~/scLMM/LMM-scRNAseq/Data/'
datafile <- paste0(dirData, "/Human_Kidney_data.rds")
data = readRDS(file = datafile)

coldata = data@meta.data
counts = GetAssayData(data, assay = 'RNA', layer = 'counts')
counts[1000:1010,1000:1010]
gene_info = data.frame(sumCount=rowSums(counts),varCount=apply(counts,1,var))
gene_info$genes = rownames(gene_info)

df_inconsist_merged = merge(df_inconsist, gene_info, by.x='genes', by.y='genes')

ggplot(df_inconsist_merged, aes(x=sclmm_tval, y=nebula_tval, color=sumCount))+
  geom_point(alpha=0.8)+theme_classic()+scale_color_viridis_b(direction = -1)
ggplot(df_inconsist_merged, aes(x=sclmm_tval, y=nebula_tval, color=varCount))+
  geom_point(alpha=0.8)+theme_classic()+scale_color_viridis_b(direction = -1)
cor(df_inconsist_merged$score, df_inconsist_merged$varCount)
cor(df_inconsist_merged$score, df_inconsist_merged$sumCount)

df_inconsist_ord


###########################################################
#################### Pathway analysis #################### 
###########################################################
source('~/RatLiver/Codes/Functions.R')
Initialize()
library(gprofiler2)
library(ggplot2)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
head(df_inconsist_ord)
num_genes = 200

enrich_res = get_gprofiler_enrich(markers=df_inconsist_ord$genes[1:num_genes], model_animal_name)
head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
title = 'inconsistent genes nebula vs scLMM'
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))


num_genes = 300
enrich_res = get_gprofiler_enrich(markers=a_cell_type_sex_df_female$genes[1:num_genes], model_animal_name)
head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
enrich_res_pos = enrich_res_pos[!is.na(enrich_res_pos$log_p),]
title = gsub(pattern = 'Male', 'Female', a_cell_type_sex)
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))




################
colnames(test)
a_cell_type_sex = 'B_cell:Male'#'Proximal_Tubule:Male' 'B_cell:Male
#a_cell_type_sex = 'CCD-like:Male'#'Proximal_Tubule:Male'
a_cell_type_sex_df = data.frame(genes= rownames(test),test[,grep(a_cell_type_sex, colnames(test))])
head(a_cell_type_sex_df)
a_cell_type_sex = gsub(':', '.', a_cell_type_sex)
a_cell_type_sex = gsub('-', '.', a_cell_type_sex)
a_cell_type_sex_df$score = -log10(a_cell_type_sex_df[[paste0(a_cell_type_sex, '_pvalue')]])* abs(a_cell_type_sex_df[[paste0(a_cell_type_sex, '_t')]])
a_cell_type_sex_df = a_cell_type_sex_df[order(a_cell_type_sex_df$score, decreasing = T),]
head(a_cell_type_sex_df,30)
a_cell_type_sex_df_male = a_cell_type_sex_df[a_cell_type_sex_df[[paste0(a_cell_type_sex, '_t')]]>0,]
head(a_cell_type_sex_df_male, 20)

num_genes = 200
enrich_res_scLMM = get_gprofiler_enrich(markers=a_cell_type_sex_df_male$genes[1:num_genes], model_animal_name)
head(enrich_res_scLMM$result,30)
enrich_res_pos_scLMM = data.frame(enrich_res_scLMM$result)
enrich_res_pos_scLMM = enrich_res_pos_scLMM[1:20,]
enrich_res_pos_scLMM = enrich_res_pos_scLMM[,colnames(enrich_res_pos_scLMM) %in% c('term_name', 'p_value')]
enrich_res_pos_scLMM$log_p = -log(enrich_res_pos_scLMM$p_value)
title = paste0(a_cell_type_sex, ' scLMM')
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))

dim(enrich_res_pos_scLMM)





a_cell_type_sex_df_nb = data.frame(genes= negbn$summary$gene,negbn$summary[,grep(a_cell_type_sex, colnames(negbn$summary))])
head(a_cell_type_sex_df_nb)
a_cell_type_sex = gsub(':', '.', a_cell_type_sex)
a_cell_type_sex = gsub('-', '.', a_cell_type_sex)
a_cell_type_sex_df_nb$score = -log10(a_cell_type_sex_df_nb[[paste0('p_',a_cell_type_sex)]])* abs(a_cell_type_sex_df_nb[[paste0('logFC_',a_cell_type_sex)]])
a_cell_type_sex_df_nb = a_cell_type_sex_df_nb[order(a_cell_type_sex_df_nb$score, decreasing = T),]
head(a_cell_type_sex_df_nb,30)
a_cell_type_sex_df_nb_male = a_cell_type_sex_df_nb[a_cell_type_sex_df_nb[[paste0('logFC_',a_cell_type_sex)]]>0,]
head(a_cell_type_sex_df_nb_male, 20)
sim(a_cell_type_sex_df_nb_male)

num_genes = 200
enrich_res_NB = get_gprofiler_enrich(markers=a_cell_type_sex_df_nb_male$genes[1:num_genes], model_animal_name)
head(enrich_res_NB$result,30)
enrich_res_pos_NB = data.frame(enrich_res$result)
enrich_res_pos_NB = enrich_res_pos_NB[1:20,]
enrich_res_pos_NB = enrich_res_pos_NB[,colnames(enrich_res_pos_NB) %in% c('term_name', 'p_value')]
enrich_res_pos_NB$log_p = -log(enrich_res_pos_NB$p_value)
title = paste0(a_cell_type_sex, ' nebula')
ggplot(enrich_res_pos_NB, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))

dim(enrich_res_NB$result)


dim(enrich_res_pos_NB)
dim(enrich_res_pos_scLMM)

sum(enrich_res_pos_NB$term_id %in% enrich_res_pos_scLMM$term_id)



