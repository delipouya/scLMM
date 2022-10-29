##Random effects by groups/clusters
##For group i:
##g(\mu_i) = X_i\beta + zi1*bi1 + ... + zik*bik = zi*bi, bik ~ N(0, sk^2)
##For all m groups (e.g., m=2):
##g(\mu) = X\beta + [z11*b11, z21*b21]^T + ... + [z1k*b1k, z2k*b2k]^T
##       = X\beta + Z1*b1 + ... + Zk*bk,  Zk = diag(z1k, z2k), bk = [b1k, b2k]^T ~ N(02, sk^2*I2)
##       = X\beta + Z*b

##designMatrix_rat_set1.rds
##inputdata_rat_set1.rds


###############################################
######### Reading the input files ###########
###############################################
##design matrix
datnm <- paste0("input_data_designMat/designMatrix_rat_set1.rds")
dat <- readRDS(file = datnm)

str(dat)
dim(dat) #[1] 23036    11
head(dat)

table(dat[, "strainLew"])
#    0     1 
#13735  9301 

all(rowSums(dat[, -(1:2)]) == 1) #[1] FALSE
table(rowSums(dat[, -(1:2)])) 
#    0     1 
# 5486 17550 
table(rowSums(dat[, -2])) 
#    1     2 
# 5486 17550 

## use clusters as fixed effects and strain as a random effect just as an example
x <- dat[, -2] #intercept and clusters
colnames(x) <- gsub("clustercluster_", "cluster", gsub("\\)", "", gsub("\\(", "", colnames(x))))
head(x)
all(x[,1] == 1)
table(rowSums(x[,-1]))
#    0     1 
# 5486 17550 
table(rowSums(x))
#    1     2 
# 5486 17550 
x <- as.matrix(x)	

z <- as.matrix(dat$strainLew) #strain as a random effect
head(z)
table(z)
#    0     1 
#13735  9301 


##raw counts
datnm <- paste0("input_data_designMat/inputdata_rat_set1_countData.rds")
dat <- readRDS(file = datnm)
dim(dat) #[1] 32883 23036
dat[1:6, 1:3]
ncells <- rowSums(dat > 0)
hist(ncells)
hist(log2(ncells))
sum(ncells >= 2^4) #[1] 12437
sum(ncells >= 20) #[1] 12103
hist(log2(ncells[ncells >= 20]))

y <- t(dat[rowSums(dat > 0) >= 20, ])
Y <- as.matrix(y)


###############################################
######### Defining the input paarameters ###########
###############################################

#family <- "binomial"
family = "poisson"
X <- x
G <- matrix(1, nrow = nrow(x), ncol = 1)
colnames(G) <- paste0("G", 1:ncol(G))
m <- 2
groupsize <- c(sum(z==0), sum(z!=0))
group <- z[,1]
group[z==0] <- 2
table(group)

##design matrix for MINQE
Z <- matrix(0, nrow = nrow(G), ncol = length(groupsize)*ncol(G))
dim(Z)
colnames(Z) <- paste0("Z", 1:ncol(Z))
m0 <- 0
for (i in 1:length(groupsize)){
  Z[(m0+1):(m0+groupsize[i]), i+(0:(ncol(G)-1))*length(groupsize)] <- G[(m0+1):(m0+groupsize[i]), ]
  m0 <- m0+groupsize[i]
}
d <- rep(length(groupsize), ncol(G))
d
head(Z)
tail(Z)


##formulas for glmer
if (family == "binomial"){
  md <- as.formula(paste0("cbind(y, size - y) ~ 0 + ", 
                          paste0(colnames(X), collapse = "+"), " + (0 + ",
                          paste0(colnames(G), collapse = "+"), " || group)"))
} else {
  md <- as.formula(paste0("y ~ 0 + ", 
                          paste0(colnames(X), collapse = "+"), " + (0 + ",
                          paste0(colnames(G), collapse = "+"), " || group)"))
}


##########
ngenes <- 20
igrp <- 0


genelist <- sample.int(ncol(Y), ngenes)
igrp <- igrp + 1

bminqe <- matrix(NA, nrow = ngenes, ncol = ncol(X))
bglmer <- bminqe
sminqe <- matrix(NA, nrow = ngenes, ncol = length(d))
sglmer <- sminqe

t1 <- Sys.time()
#for (i in 1:ngenes){
for (i in 1:length(genelist)){
  y <- Y[, genelist[i]]
  
  ##lme4::glmer
  dat <- data.frame(y, X, G, group = group)
  fitglmer <- glmer(md, data = dat, family = family)
  #summary(fitglmer)
  #v <- vcov.merMod(fitglmer, corr=TRUE)
  #as(v, "corMatrix")
  #fixef(fitglmer)
  #as.data.frame(VarCorr(fitglmer))[, "vcov"]
  
  #est <- glmm.minqe(y, X, Z, d = d, family = family, phi = phi, method = "I", plot.it = T)
  est <- glmm.minqe(y, X, Z, d = d, family = family, niter = 10, niter.minqe = 5, nsamples = 100, phi = phi, method = "I", plot.it = F)
  
  bglmer[i, ] <- fixef(fitglmer)
  sglmer[i, ] <- as.data.frame(VarCorr(fitglmer))[, "vcov"]
  bminqe[i, ] <- est$beta
  sminqe[i, ] <- est$theta[1:length(d)]
  
}##FALSE
t2 <- Sys.time()
difftime(t2, t1)
##