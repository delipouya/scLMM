""" Dimension reduction of complex data with supervision from auxiliary information. 
The package contains a series of methods for different data types (e.g., multi-view or multi-way data) 
including the supervised singular value decomposition (SupSVD), supervised sparse and functional principal component (SupSFPC), 
supervised integrated factor analysis (SIFA) and supervised PARAFAC/CANDECOMP factorization (SupCP). When auxiliary data 
are available and potentially affect the intrinsic structure of the data of interest, the methods will accurately recover 
the underlying low-rank structure by taking into account the supervision from the auxiliary data.
"""

#install.packages("https://cran.r-project.org/src/contrib/Archive/SuperPCA/SuperPCA_0.4.0.tar.gz", repos=NULL)
#BiocManager::install('matlabr') #‘fBasics’, ‘R.matlab’, ‘matlabr’ 

## source code: https://rdrr.io/cran/SuperPCA/f/
## API and list of functions: https://rdrr.io/cran/SuperPCA/api/

library(SuperPCA)
help('SuperPCA')


###### Supervised SVD tutorial
r=2
Y <- matrix(rnorm(400,0,1),nrow=100) 
dim(Y)
head(Y)
colnames(Y) = paste0('cov_', 1:ncol(Y))
rownames(Y) = paste0('sample_', 1:nrow(Y))

B <- c(-1,1,-sqrt(3/2),-1)
B <- cbind(B,c(1,1,-1,sqrt(3/2)))
V <- matrix(rnorm(68*2),68,2)
dim(V)
Fmatrix <- matrix(MASS::mvrnorm(n=1*100,rep(0,2),matrix(c(9,0,0,4),2,2)),100,2)
dim(Fmatrix)
E <- matrix(rnorm(100*68,0,3),100,68)
Yc <- scale(Y,center=TRUE,scale=FALSE)

# Case 1 (supsvd) X = YBV^T+FV^T+E
X1 <- Y%*%tcrossprod(B,V)+tcrossprod(Fmatrix,V)+E
X1c <- scale(X1,center=TRUE,scale=FALSE)
res1 = SupPCA(Yc,X1c,r)
names(res1)


#  Case 2 (PCA) X = FV^T+E
X2 <- tcrossprod(Fmatrix,V)+E
X2c <-scale(X2,center=TRUE,scale=FALSE)
res2 = SupPCA(Yc,X2c,r)
names(res2)
# Case 3 (RRR) X = YBV^T+E
X3 <- Y%*%tcrossprod(B,V)+E
X3c <- scale(X3,center=TRUE,scale=FALSE)
res3 = SupPCA(Yc,X3c,r)
names(res3)
res3$U

###### SIFA tutorial
r0 <- 2
r <- c(3,3)
V <- matrix(stats::rnorm(10*2),10,2)
Fmatrix <- matrix(MASS::mvrnorm(n=1*500,rep(0,2),matrix(c(9,0,0,4),2,2)),500,2)
E <- matrix(stats::rnorm(500*10,0,3),500,10)
X <- tcrossprod(Fmatrix,V)+E
X <-scale(X,center=TRUE,scale=FALSE)
Y1 <- matrix(stats::rnorm(500*200,0,1),500,200)
Y2 <- matrix(stats::rnorm(500*200,0,1),500,200)
Y <- list(Y1,Y2)
SIFA(X,Y,r0,r,max_niter=200)


###### SupSFPCA tutorial
library(spls)
data(yeast)
r <- 4
ydata <- as.data.frame(yeast[1])
xdata <- as.data.frame(yeast[2])
yc <- scale(ydata,center = TRUE,scale=FALSE)
xc <- scale(xdata,center=TRUE,scale=FALSE)
SupSFPCA_res = SupSFPCA(yc,xc,r)
names(SupSFPCA_res)

