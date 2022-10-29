source('~/RatLiver/Codes/Functions.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(lme4) # load library
library(arm) # convenience functions for regression in R
library(dplyr)

gene_exp_biased = read.table('~/scLMM/biased_model_inputs/Y_div.txt')
gene_exp_biased = read.table('~/scLMM/biased_model_inputs/Y_div_2.txt')
Z = read.table('~/scLMM/biased_model_inputs/Z_biased.txt')
head(Z)
head(gene_exp_biased)
Z = Z / ncol(Z)


random_effect_chr = ifelse(Z$V1==1, 'A', 'B')
table(random_effect_chr)
table(Z$V1)

input.data = data.frame(y=gene_exp_biased$V1, 
                        kernel=as.factor(random_effect_chr)
)
head(input.data)
attach(input.data)

par(mfrow=c(1,1))
barplot(table(input.data$kernel), ylab = "Frequency", main = "kernel")
ggplot(input.data, aes(x=kernel, y=y))+geom_boxplot()+theme_classic()


##############################################
Model = lmer(y ~ (1|kernel), data=input.data, REML=T) 
getME(Model, "X")
getME(Model, "Z")

summary(Model)
df_3 = getResDF(Model_3)
ggplot(df_3, aes(x=fitted, y=residuals, color=strain))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()
ggplot(df_3, aes(x=fitted, y=residuals, color=cluster))+
  geom_point(alpha=0.3, size=2.5)+theme_classic()+
  scale_color_manual(values = c('black', 'red3'))
hist(df_3$residuals, main='residuals',xlab='residuals')
qqnorm(residuals(Model_3))
coef(Model_3)
