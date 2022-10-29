lmer(gene~ cluster +(1|strain), data=input.data, REML=F) 
a=rep(1:2,3)
b=rep(1:3,2)
x=1:6
df=data.frame(A=a,B=b,x=x)
# Lie and pretend there's a level 0 in each factor.
df$A=factor(a,as.character(0:2))

df$B=factor(b,as.character(0:3))

mm=model.matrix (~B + (1|A) ,df) #B:x

print(mm)



sig_1 = 0.2
sig_2 = 0.5
sig_res = 0.3
num_samples = 1000
num_feat = 10
num_simulate = 50

title = paste0('num-samples: ', num_samples, ', num-features ', num_feat, ', num-simulations: ', num_simulate)
subtitle =  paste0('sig1: ', sig_1, ', sig2: ', sig_2, ', sig-residual: ', sig_res)

df = read.csv('~/scLMM/gene1_simulation_1000_numFeat10.csv')
df.m = reshape2::melt(df[,1:3])

files = list.files(pattern = 'gene_*', '~/scLMM/python_codes/sim_50_numFeat_10_numSamp_1000/',full.names = T)
files = list.files(pattern = 'gene_*', '~/scLMM/python_codes/sim_50_numFeat_10_numSamp_10000/',full.names = T)

inputs = lapply(files, read.csv)
names(inputs) = paste0('gene_', 1:length(files)) 
lapply(inputs, function(df) table(df$convergence))

inputs = sapply(1:length(inputs), 
                function(i) {df = data.frame(inputs[[i]]); df$feature=names(inputs)[i];  reshape2::melt(df[,c(1:3,5)])},
                simplify = F)



merged =  do.call(rbind,inputs)
head(merged)

ggplot(merged, aes(x=feature, y=value, fill=variable))+geom_boxplot()+
  geom_hline(yintercept=sig_1, linetype="dashed", color = "red")+
  geom_hline(yintercept=sig_2, linetype="dashed", color = "green")+
  geom_hline(yintercept=sig_res, linetype="dashed", color = "blue")+
  theme_bw()+ylab('sigma')+labs(title = title, subtitle = subtitle)

ggplot(merged, aes(x=variable, y=value, fill=feature))+geom_boxplot()+
  geom_hline(yintercept=sig_1, linetype="dashed", color = "red")+
  geom_hline(yintercept=sig_2, linetype="dashed", color = "green")+
  geom_hline(yintercept=sig_res, linetype="dashed", color = "blue")+
  theme_bw()+ylab('sigma')+labs(title = title, subtitle = subtitle)



df1 = merged[merged$variable =='sig1',]
ggplot(df1, aes(x=feature, y=value, fill=feature))+geom_boxplot()+
  geom_hline(yintercept=sig_1, linetype="dashed", color = "blue")+theme_bw()+ylab('sigma-1')

df2 = merged[merged$variable =='sig2',]
ggplot(df2, aes(x=feature, y=value, fill=feature))+geom_boxplot()+
  geom_hline(yintercept=sig_2, linetype="dashed", color = "blue")+theme_bw()+ylab('sigma-2')

df3 = merged[merged$variable =='sig_res',]
ggplot(df3, aes(x=feature, y=value, fill=feature))+geom_boxplot()+
  geom_hline(yintercept=sig_res, linetype="dashed", color = "blue")+theme_bw()+ylab('sig_residuals')





