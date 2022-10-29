source('~/RatLiver/Codes/Functions.R')
source('~/RatLiverCodes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(lme4) # load library
library(arm) # convenience functions for regression in R
library(dplyr)

############
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)

sample_size = 1000
sample_size = ncol(merged_samples)
selected_umi_index = sample(size = sample_size, 1:ncol(merged_samples))

feature_size = 500
feature_size = nrow(merged_samples)
selected_gene_index = sample(size = feature_size, 1:nrow(merged_samples))

#Itgal_index = which(rownames(merged_samples) == 'Itgal')
#data_sub = GetAssayData(merged_samples[c(selected_gene_index, Itgal_index),selected_umi_index])
merged_samples_org = merged_samples
merged_samples = merged_samples[,merged_samples$cluster %in% c(5, 10)]
data_sub = GetAssayData(merged_samples)
input.data = data.frame(gene=data_sub['Itgal',], #'Itgal'
                        cluster=as.factor(merged_samples$cluster), #[selected_umi_index]
                        sample=as.factor(merged_samples$sample_name),
                        strain=as.factor(merged_samples$strain),
                        numFeatures=merged_samples$nFeature_RNA,
                        librarySize=merged_samples$nCount_RNA
                        )

attach(input.data)
input.data$gene = log10(1 + input.data$gene_2)
head(input.data)
dim(input.data)
hist(input.data$gene, main ='Itgal',xlab='Expression values', breaks = seq(from = 0,to=2,by=.2))

par(mfrow=c(1,1))
barplot(table(input.data$cluster), ylab = "Frequency", main = "Cluster")
barplot(table(input.data$strain), ylab = "Frequency", main = "Strain")

ggplot(input.data, aes(x=strain, y=gene))+geom_boxplot()+theme_classic()
ggplot(input.data, aes(x=cluster, y=gene))+geom_boxplot()+theme_classic()
ggplot(input.data, aes(x=cluster, y=log10(gene+1)))+geom_boxplot()+theme_classic()



getResDF <- function(Model){
  df = data.frame(fitted=fitted(Model),
                  residuals=residuals(Model),
                  sample=merged_samples$sample_name,
                  strain=merged_samples$strain,
                  cluster=merged_samples$cluster)
  return(df)
}


##############################################
Model_0 = lm(gene~ strain, data=input.data) 
summary(Model_0)
df_0 = getResDF(Model_0)
ggplot(df_0, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_0, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_0$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_0))
table(df$fitted, df$residuals)

##############################################
Model_1 = lm(gene ~ cluster, data=input.data) 
summary(Model_1)
df_1 = getResDF(Model_1)
ggplot(df_1, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_1, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_1$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_1))

##############################################
Model_2 = lm(gene~ strain+cluster, data=input.data) 
summary(Model_2)
df_2 = getResDF(Model_2)
ggplot(df_2, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_2, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_2$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_2))

anova(Model_0, Model_1, Model_2)
#Model_2g = gls(gene~ strain+cluster, data=input.data, method='ML') 

##############################################
Model_3 = lmer(gene~ strain +(1|cluster), data=input.data, REML=F) 
getME(Model_3, "X")
getME(Model_3, "Z")

summary(Model_3)
df_3 = getResDF(Model_3)
ggplot(df_3, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_3, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_3$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_3))
coef(Model_3)

##############################################
Model_4 = lmer(gene~ cluster +(1|strain), data=input.data, REML=F) 
summary(Model_4)
df_4 = getResDF(Model_4)
ggplot(df_4, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_4, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_4$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_4))
coef(Model_4)

##############################################
### modeling by-cluster variablity in how strain effects gene expression 
Model_5 = lmer(gene~ strain +(1+strain|cluster), data=input.data, REML=F) 
summary(Model_5)
df_5 = getResDF(Model_5)
ggplot(df_5, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_5, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_5$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_5))
coef(Model_5)

##############################################
Model_allRand = lmer(gene~ 1+(1|strain)+(1|cluster), data=input.data, REML=F) 
summary(Model_allRand)
df_rand = getResDF(Model_allRand)
ggplot(df_rand, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_rand, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_rand$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_allRand))
coef(Model_allRand)

##############################################
Model_6 = lmer(gene~ strain +(1+strain|cluster)+(1|sample), data=input.data, REML=F) 
df_6 = getResDF(Model_allRand)
ggplot(df_6, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_6, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_6$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_6))
coef(Model_6)

##############################################
Model_7 = lmer(gene~ cluster + (1|strain/sample), data=input.data, REML=F)
df_7 = getResDF(Model_7)
ggplot(df_7, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_7, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_7$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_7))
coef(Model_7)

##############################################
Model_8 = lmer(gene~ strain + (1+strain|cluster) + (1|strain/sample), data=input.data, REML=F)
summary(Model_8)
df_8 = getResDF(Model_8)
ggplot(df_8, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_8, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_8$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_8))
coef(Model_8)

##############################################
Model_9 = lmer(gene~ strain + cluster + (1|sample), data=input.data, REML=F)
summary(Model_9)
df_9 = getResDF(Model_9)
ggplot(df_8, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_8, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_8$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_9))
coef(Model_9)

anova(Model_3, Model_4, Model_5, Model_allRand, Model_6, Model_7, Model_8, Model_9)


plot(input.data$gene, df_8$fitted, main='Model 8')
plot(input.data$gene, df_7$fitted, main='Model 7')
plot(input.data$gene, df_6$fitted, main='Model 6')
plot(input.data$gene, df_5$fitted, main='Model 5')
par(mfrow=c(2,2))
plot(input.data$gene, df_4$fitted, main='Model 4')
plot(input.data$gene, df_3$fitted, main='Model 3')
plot(input.data$gene, df_2$fitted, main='Model 2')
plot(input.data$gene, df_1$fitted, main='Model 1')
par(mfrow=c(2,2))


##########################################

gene_indices = match(c('Itgal', 'Cth', 'Rgs1', 'Cyp4a2'),
                     rownames(merged_samples))
input.data = data.frame(t(data_sub[c(1:100,gene_indices),]),
                        cluster=as.factor(merged_samples$cluster), #[selected_umi_index]
                        sample=as.factor(merged_samples$sample_name),
                        strain=as.factor(merged_samples$strain),
                        numFeatures=merged_samples$nFeature_RNA,
                        librarySize=merged_samples$nCount_RNA)

input.data = data.frame(t(data_sub),
                        cluster=as.factor(merged_samples$cluster), #[selected_umi_index]
                        sample=as.factor(merged_samples$sample_name),
                        strain=as.factor(merged_samples$strain),
                        numFeatures=merged_samples$nFeature_RNA,
                        librarySize=merged_samples$nCount_RNA)
head(input.data)
dim(input.data)
dim(t(data_sub))

model = lm( as.matrix(input.data[,c(1:104)]) ~ strain + cluster, data=input.data)
model = lm( as.matrix(input.data[,c(1:12004)]) ~ strain + cluster, data=input.data) 
summary(model)
model.coef = data.frame(t(coef(model)))
dim(model.coef)
head(model.coef)

head(model.coef[order(model.coef$strainLew, decreasing = T),],15)
head(model.coef[order(model.coef$strainLew, decreasing = F),],15)

head(model.coef[order(model.coef$cluster5, decreasing = T),], 15)
head(model.coef[order(model.coef$cluster5, decreasing = F),], 15)


# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
# https://rpubs.com/bbolker/3336
# https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/
# https://medium.com/analytics-vidhya/performing-multivariate-mixed-modeling-7c925f015f39
library(MCMCglmm)
m1_idyr <- MCMCglmm(cbind(gene_1, gene_2) ~ strain + (1|cluster),
                    data = input.data, family = c("gaussian", "gaussian"))

MCMCglmm(cbind(y.hol, y.car) ~ trait - 1, random = ~us(trait):id, 
         rcov = ~us(trait):units, data = Spending, family = c("gaussian", "gaussian"), 
         verbose = FALSE)



MCMCglmm(cbind(Ulbp1,Lrp11)  ~ strain, random = ~cluster, data = input.data, verbose = FALSE)











##########################################
res <- model.matrix(~cluster, data = input.data)
head(res[, -1])
contrasts(input.data$cluster)
contrasts(input.data$strain)

library(jtools)
summ(test_m)
ranova(M1)

library(reghelper)
simple_slopes(test_m)
graph_model(test_m,y=strain,x=cluster,lines=sample)

display(MLexamp)
AIC(MLexamp)

# Install (if necessary) and load nlme and lme4
library(nlme)
library(lme4)
# Load dataset, inspect size and additional info
data(Arabidopsis)
dim(Arabidopsis) # 625 observations, 8 variables
attach(Arabidopsis)

#### fitting a GML (does not include random effects) -> proper null model with respect to random effects. 
GLM <- gls(total.fruits ~ rack + nutrient + amd + status,
           method = "ML")
summary(GLM)
#######################################
#Let’s check how the random intercepts and slopes distribute in the highest level (i.e. gen within popu).
plot(ranef(lmm6.2, level = 2))


