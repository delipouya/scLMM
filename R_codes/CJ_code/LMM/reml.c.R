##################################################
##Restricted maximum likelihood (REML)  for LMM estimation 
##Iterative algorithms for REML: 
##(1) Gradient methods: Newton–Raphson (NR), Fisher scoring (FS), and average information (AI)
##(2) Expectation-maximization (EM) algorithm, and 
##(3) Iterated MINQE (minimum norm quadratic estimation) 
##
##reml.c: 
##- REML algorithm working on columns of data matrix
##- Assuming the number of columns is not large
##- REML with FS
##################################################

reml.c <- function(Y, X, Z, d, s0 = NULL, max.iter = 50, epsilon = 1e-3) {
##Inputs
##Y: a matrix of observed measurements
##X:  a design matrix for fixed effects 
##Z = [Z1, ..., Zk],  a design matrix for random effects 
##  k: number of regions (kernels),
##d = (m1,...,mk), 
##  mi = ncol(Zi), number of columns in Zi
##  m1 + ... + mk = ncol(Z), number of columns in Z	
##s = (s1, ...,sk, s_{k+1}), a vector of variance components
##  si = sigma_i^2, i = 1, ..., k
##  s_{k+1} = sigma_0^2
##  s0, initial values of the variance components
##
##Outputs
##theta: a matrix of the variance component estimates for each sample
##coef: a matrix of the fixed effects
##t: t-values of the fixed effects
##pvalue: p-values of the fixed effects
##Wald.stat: Wald statistics of the fixed effects
##Wald.pvalue: Wald test p-values for the fixed effects


Y <- as.matrix(Y)
X <- as.matrix(X)
Z <- as.matrix(Z)
n <- nrow(Y) 
p <- ncol(X)
k <- length(d)  

stopifnot(sum(d) == ncol(Z))

##########
##estimated variance components in iterations
##dl: first derivatives of log likelihood
##xx = (X'X)^{-1}
zz <- t(Z)%*%Z
yy <- colSums(Y*Y) 
xx <- ginv(t(X)%*%X)
xy <- t(X)%*%Y
X <- t(Z)%*%X
xxz <- xx%*%t(X)
Y <- t(Z)%*%Y

##zrz = Z'RZ
##zry = Z'Ry
##yry = [y1'Ry1,...,ym'Rym]
zrz <- zz - X%*%(xx%*%t(X))
zry <- Y - X%*%(xx%*%xy)
yry <- yy - colSums(xy*(xx%*%xy))


##########
niter <- NULL
theta <- NULL
pvtheta <- NULL
stat <- NULL
pvstat <- NULL
beta <- NULL
pvbeta <- NULL
tvbeta <- NULL
dlogL <- NULL

for (jy in 1:ncol(Y)) {
	
	if (is.null(s0)) {
		s <- yry[jy]/(n-p)
		s <- c(rep(0, k), s)
	} else s <- s0

dl <- 100
iter <- 0

while ((max(abs(dl)) > epsilon)	& (iter < max.iter)){
	iter <- iter + 1
	
	fs <- matrix(NA, k+1, k+1)	##Fisher scoring matrix
	dl <- rep(NA, k+1)
		
	sr <- s[1:k]/s[k+1]
	##M = (SZ'RZ + I)^{-1}
	M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
	ZRZ <- zrz%*%M
	ZR2Z <- ZRZ%*%M
	yRZ <- t(zry[, jy])%*%M
	
	mi <- 0
	for (i in 1:k){	
		ik <- (mi+1):(mi+d[i])
		dl[i] <- (sum((yRZ[ik])^2)/s[k+1]^2 - sum(diag(ZRZ[ik, ik]))/s[k+1])/2

	mj <- 0
	for (j in 1:i){
		ji <- (mj+1):(mj+d[j])
		fs[i, j] <- sum((ZRZ[ji, ik])^2)/s[k+1]^2/2
		fs[j, i] <- fs[i, j]
		mj <- mj + d[j]
		}
		
	j <- k+1		
	fs[i, j] <- sum(diag(ZR2Z[ik, ik]))/s[k+1]^2/2
	fs[j, i] <- fs[i, j]
	mi <- mi + d[i]
	}
	
	i <- k+1
	fs[i, i] <- (n - p - sum(d) + sum(t(M)*M))/s[k+1]^2/2
	
	yR2y <- yry[jy] - sum(((t(M) + diag(sum(d)))%*%zry[, jy])*(M%*%(rep(sr, times = d)*zry[, jy])))
	dl[i] <-  (yR2y/s[k+1]^2 - (n-p-sum(d)+sum(diag(M)))/s[k+1])/2

	##The H, FS, and AI matrices can be singular.
	#Minv <- solve(M)
	Minv <- ginv(fs)
	
	s <- s + Minv%*%dl
	#s <- s1*(s1 >= 0) + s*(s1 < 0)
	#s <- s*(s>=0)
	#theta.iter <- cbind(theta.iter, s)
	#dl.iter <- cbind(dl.iter, dl)
	}
	
##
if (max(abs(dl)) > epsilon) warning(paste0("One of the first derivatives of log likelihood for the column ", jy, " > epsilon, ", epsilon))


##########
##Wald test: beta = 0
##b: beta = (X'V1X)^{-1}X'V1y
##Vb: var(beta) = (X'V1X)^{-1}
##W = beta'*Var(beta)^{-1}*beta ~ chisq(p), asymptotically
##W = y'V1X(X'V1X)^{-1}X'V1y, V1 = V^{-1}.
##v1x = V^{-1}X = V1X

sr <- s[1:k]/s[k+1]
M <- solve(sweep(zz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*") 
xvx <- xx + xxz%*%(ginv(diag(sum(d)) - M%*%(X%*%xxz))%*%(M%*%t(xxz)))
xvy <- xy[, jy] - t(X)%*%(M%*%Y[, jy])
b <- xvx%*%xvy
tv <- b/sqrt(abs(diag(xvx)))/s[k+1]
pvb <- 2*pt(-abs(tv), df = n - p)
W <- sum(xvy*b)/s[k+1]
pvw <- 1 - pchisq(W, df = p)

##outputs
niter <- c(niter, iter)
theta <- cbind(theta, s)
#pvtheta <- cbind(pvtheta, 1 - pchisq(n*abs(s)/diag(solve(fs)), df = 1))
stat <- c(stat, W)
pvstat <- c(pvstat, pvw)
beta <- cbind(beta, b)
tvbeta <- cbind(tvbeta, tv)
pvbeta <- cbind(pvbeta, pvb)
dlogL <- cbind(dlogL, dl)
}


list(method = "FS", dlogL = dlogL, niter = niter, coef = beta, t = tvbeta, 
     pvalue = pvbeta, Wald.stat = stat, Wald.pvalue = pvstat, theta = theta)
}

##################################################
##source("R/reml.c.R")
##a <- reml.c(Y, X, Z, d = c(2,ncol(Z)-2))
##b <- reml(Y, X, Z, d = c(2,ncol(Z)-2))
##
##a <- reml.c(Y, X, Z, d = ncol(Z))
##a <- reml.c(Y, X, Z, d = ncol(Z), max.iter = 10)
##b <- reml(Y, X, Z, d = ncol(Z), niter = 10)
##  range(a$FS-b$FS)
##  range(a$theta.iter-b$theta.iter)
##  range(a$dl.iter-b$dl.iter)
##################################################
