##################################################
##Restricted maximum likelihood (REML)  for LMM estimation 
##Iterative algorithms for REML: 
##(1) Gradient methods: Newton–Raphson (NR), Fisher scoring (FS), and average information (AI)
##(2) Expectation-maximization (EM) algorithm, and 
##(3) Iterated MINQE (minimum norm quadratic estimation) 
##
##reml.c: (v5.2)
##- REML algorithm working on columns of data matrix
##- Assuming the number of columns is not large
##- REML with FS
##
##lmmfitSS: using correlation summary statistics (SS), XY, ZX and ZY, to fit LMM.
##
##Inputs
##XXinv = (X'X)^{-1}
##ZZ: ZZ <- t(Z)%*%Z
##Ynorm: Ynorm <- colSums(Y*Y) 
##XXinv: XXinv <- ginv(t(X)%*%X)
##XY: XY <- t(X)%*%Y
##ZX: ZX <- t(Z)%*%X
##ZY: ZY <- t(Z)%*%Y
#n <- nrow(Y) 
##
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
##################################################

lmmfitSS <- function(XY, ZX, ZY, ZZ, XXinv, Ynorm, n, d, s0 = NULL, 
		method = "REML-FS", max.iter = 50, epsilon = 1e-5)
{
#XY <- as.matrix(XY)
#ZX <- as.matrix(ZX)
#ZY <- as.matrix(ZY)
#ZZ <- as.matrix(ZZ)
#XXinv <- as.matrix(XXinv)

p <- ncol(ZX)
k <- length(d)  

stopifnot(sum(d) == ncol(ZZ))

##zrz = Z'RZ
##zry = Z'Ry
##yry = [y1'Ry1,...,ym'Rym]
xxz <- XXinv%*%t(ZX)
zrz <- ZZ - ZX%*%(XXinv%*%t(ZX))
zry <- ZY - ZX%*%(XXinv%*%XY)
yry <- Ynorm - colSums(XY*(XXinv%*%XY))


##########
##estimated variance components in iterations
##dl: first derivatives of log likelihood

niter <- NULL
dlogL <- NULL
theta <- matrix(nrow = k + 1, ncol = ncol(XY), dimname = list(c(1:k, 0), colnames(XY)))
beta <- matrix(nrow = nrow(XY), ncol = ncol(XY), dimnames = dimnames(XY))
covbeta <- array(dim = c(nrow(XY), nrow(XY), ncol(XY)), 
	dimnames = list(rownames(XY), rownames(XY), colnames(XY)))

for (jy in 1:ncol(ZY)) {
	
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
		dl[i] <- (sum((yRZ[ik])^2)/s[k+1]^2 - sum(diag(ZRZ[ik, ik, drop = FALSE]))/s[k+1])/2

	mj <- 0
	for (j in 1:i){
		ji <- (mj+1):(mj+d[j])
		fs[i, j] <- sum((ZRZ[ji, ik])^2)/s[k+1]^2/2
		fs[j, i] <- fs[i, j]
		mj <- mj + d[j]
		}
		
	j <- k+1		
	fs[i, j] <- sum(diag(ZR2Z[ik, ik, drop = FALSE]))/s[k+1]^2/2
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
	}
	
##
if (max(abs(dl)) > epsilon) warning(paste0("One of the first derivatives of log likelihood for the column ", jy, " > epsilon, ", epsilon))


##########
##Wald test: beta = 0
##b: beta = (X'V1X)^{-1}X'V1y
##Vb: var(beta) = (X'V1X)^{-1}
##W = beta'*Var(beta)^{-1}*beta ~ chisq(p), asymptotically
##W = y'V1X(X'V1X)^{-1}X'V1y, V1 = V^{-1}.
##v1x = V^{-1}ZX = V1X

sr <- s[1:k]/s[k+1]
M <- solve(sweep(ZZ, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*") 
xvx <- XXinv + xxz%*%(ginv(diag(sum(d)) - M%*%(ZX%*%xxz))%*%(M%*%t(xxz)))
xvy <- XY[, jy] - t(ZX)%*%(M%*%ZY[, jy])
b <- xvx%*%xvy
covbeta[,,jy] <- (xvx + t(xvx))*(s[k+1]/2)
##xvx*s[k+1]

##t-value
#tv <- b/sqrt(abs(diag(xvx)))/s[k+1]
#tv <- b/sqrt(abs(diag(xvx))*s[k+1])
#pvb <- 2*pt(-abs(tv), df = n - p)
##Wald test
#W <- sum(xvy*b)/s[k+1]
#pvw <- 1 - pchisq(W, df = p)

##z-score
#zs <- tv
#tNeg <- (tv < 0)
#zs[tNeg] <- qnorm(pt(tv[tNeg], df = n-p))
#zs[!tNeg] <- - qnorm(pt(-tv[!tNeg], df = n-p))

##outputs
niter <- c(niter, iter)
theta[, jy] <- s 
beta[, jy] <- b 
dlogL <- cbind(dlogL, dl)
}

list(method = method, dlogL = dlogL, niter = niter, coef = beta, cov = covbeta, df = n-p, theta = theta)
}
