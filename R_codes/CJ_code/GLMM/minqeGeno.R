
##################################################
##MINQE computing from simplified Rao's formulas with 
##- (1) trace estimation 
##- (2.1) conjugate gradients using R function cg.mat or
##- (2.2) Jacobi (diagonal) preconditioned conjugate gradients using R function pcgDiag 
##
##inputs
##y: a vector of observations
##X: covariates matrix
##G: genotypes matrix G = [G1, ... , G_p], Gi is a matrix of genotypes in a region.
##d: numbers of genotypes in each G_i, i = 1,...,p. 
##nsamples: number of samples for stochastic trace estimation
##epsilon: positive convergence tolerance in conjugate gradient algorithm
##
##outputs
##- eta: prediction of y
##- beta: estimated fixed effects 
##- theta: estimated covariance components, i.e, sigma squared
##- theta.iter: estimated covariance components in each iteration
##################################################

minqeGeno <- function(y, X, G, d, theta = rep(1, length(d)+1), phi = NULL, method = c("IU", "I"), niter = 1,
	nsamples = 100, sampleDist = rnorm, epsilon = 1e-8, plot.it = FALSE)
{
##d = (m1,...,mk), numbers of geneotypes
##Using Rademacher distribution (Hutchinson's estimator) has lower variance compared to the Gaussian.
##Trace estimation converges fast as total number of genotypes becomes larger.
##Hutchinson's estimator:
##Rademacher distribution = 2*Bernoulli - 1
##sampleDist <- function(n) 2*rbinom(n, size = 1, prob = 0.5) - 1

method <- match.arg(method)	
	k <- length(d)  #number of regions
	n <- length(y)
	p <- ncol(X)
	m <- sum(d)    #total number of SNPs
	
	if (!is.null(phi)) theta[k+1] <- phi
	s <- theta[1:k]/theta[k+1]
	
	##########	
	##MINQE(I,U):
	##R = I_n - X(X^TX)^{-1}X^T, tr(R) = n-p
	##MINQE(I):	
	##R = I_n, tr(R) = n
	
	G0 <- G
	
	trR <- n
	##G:=RG
	if (method == "IU"){
		G <- G - X%*%solve(t(X)%*%X)%*%(t(X)%*%G)
		trR <- n - p
	}
	
	##########
	theta.iter <- matrix(NA, nrow = niter, ncol = k+1)
	colnames(theta.iter) <- paste0("theta", 1:(k+1))
	rownames(theta.iter) <- 1:niter
	
	for (iter in 1:niter){
	
	##########	
	##Compute beta and eta
	##(A'A+phi*I)x = b
	##A = G0*D
	##phi = theta[k+1]
	A <- NULL
	mi <- 0
	for (i in 1:k){
		A <- cbind(A, sqrt(s[i])*G0[, (mi+1):(mi+d[i]), drop = FALSE])
		mi <- mi + d[i]
		}

	VX <- pcgDiag(A, b = t(A)%*%X, lambda = theta[k+1], epsilon = epsilon, verbose = F)	
	VX <- X - A%*%VX
	beta <- solve(t(X)%*%VX, t(VX)%*%y)
	
	e <- y - X%*%beta
	Ve <- pcgDiag(A, b = t(A)%*%e, lambda = theta[k+1], epsilon = epsilon, verbose = F)	
	Ve <- e - A%*%Ve
	eta <- y - Ve
	
	##########	
	##computing matrix S
	##- conjugate gradients using R function pcgDiag: 
	##  (A'A+I)x = b or x = (A'A+I)^{-1}b
	##  A = RGD
	
	A <- NULL
	mi <- 0
	for (i in 1:k){
		A <- cbind(A, sqrt(s[i])*G[, (mi+1):(mi+d[i]), drop = FALSE])
		mi <- mi + d[i]
		}
	
	##matrix consisting of random m-vectors
	z <- matrix(sampleDist(m*nsamples), nrow = m, ncol = nsamples)	

	##S
	S <- matrix(NA, k+1, k+1)
	
	##S[k+1, k+1]
	Mz <- pcgDiag(A, b = z, epsilon = epsilon, verbose = F)	
	fz <- colSums(Mz^2)
	S[k+1,k+1] <- trR - m + mean(fz)
	
	##S[i,j]
	mi <- 0
	for (i in 1:k){
	mj <- 0			
	i12 <- (mi+1):(mi+d[i])	
	b <- t(A)%*%(G[, i12, drop = FALSE] %*% z[i12, ])
	Mb <- pcgDiag(A, b = b, epsilon = epsilon, verbose = F)
	for (j in 1:i){
		j12 <- (mj+1):(mj+d[j])		
		Mbj <- t(G[, j12])%*%((G[, i12, drop = FALSE]%*%z[i12, ]) - A%*%Mb)
		fz <- colSums(Mbj^2)
		S[i,j] <- mean(fz)
		S[j,i] <- S[i,j]
		mj <- mj + d[j]
	}
	Mbi <- t(G[, i12])%*%((G[, i12, drop = FALSE]%*%z[i12, ]) - A%*%Mb)
	fz <- colSums(z[i12, ]*Mbi) - colSums(Mb^2)			
	S[k+1,i] <- mean(fz)
	S[i,k+1] <- S[k+1,i]
	mi <- mi + d[i]
	}
	
	##########
	##computing vector q	
	##- conjugate gradients: (A'A+I)x = b or x = (A'A+I)^{-1}b
	##M = ((DG'R)(RGD) + I_n)^{-1}
	##For MINQE(I), the matrix M in S is different from the matrix M in q.
	##G:=RG
	
	if (method == "I"){
		G <- G - X%*%solve(t(X)%*%X)%*%(t(X)%*%G)
		}
	
	q <- rep(NA, k+1)
	
	##gy = DG'Ry
	##A = RGD
	A <- NULL
	gy <- t(G)%*%y #G'Ry
	mi <- 0
	for (i in 1:k){
		A <- cbind(A, sqrt(s[i])*G[, (mi+1):(mi+d[i]), drop = FALSE])
		gy[(mi+1):(mi+d[i])] <- sqrt(s[i])*gy[(mi+1):(mi+d[i])]
		mi <- mi + d[i]
	}

	My <- pcgDiag(A = A, b = gy, epsilon = epsilon) #M*dgry
	mi <- 0
	for (i in 1:k){	
		i12 <- (mi+1):(mi+d[i])
		myi <- t(G[, i12, drop = FALSE])%*%(y - A%*%My)
		q[i] <- sum(myi^2)
		mi <- mi + d[i]
	}
	xy <- t(X)%*%y
	yry <- sum(y*y) - sum(xy*(solve(t(X)%*%X)%*%xy))
	q[k+1] <- yry - sum(gy*My) - sum(My^2)
	
	#theta <- c(solve(S)%*%q)
	##library(nnls)
	if (is.null(phi)) theta <- nnls(S, q)$x
	else theta[1:k] <- nnls(as.matrix(S[-(k+1), -(k+1)]), q[-(k+1)] - theta[k+1]*S[-(k+1), k+1])$x
	#else theta[1:k] <- nnls(as.matrix(S[, -(k+1)]), q - theta[k+1]*S[, k+1])$x
	
	if (theta[k+1] == 0) s <- theta[1:k]
	else s <- theta[1:k]/theta[k+1]
	
	
	##non-negative least squares
	##method 1
	##library(nnls)
	#theta <- nnls(S, q)$x
	#theta <- nnls(Sc, q)$x
	##method 2
	##library(glmnet)
	#theta <- glmnet(S, q, lambda = 0, lower.limits = 0, intercept = FALSE)$beta
	#theta <- as.vector(theta)

	theta.iter[iter, ] <- theta	
}

if (plot.it){
	plot(theta.iter[,1], ylim = range(theta.iter), xlab = "Iteration", ylab = "theta", type = "n")
	colset <- rainbow(max(7, k+1))
	for (j in 1:ncol(theta.iter)) {
		lines(theta.iter[,j], col = colset[j])
		points(theta.iter[,j], col = colset[j], pch = 16, cex = 0.8)
		}
}

list(theta.iter = theta.iter, eta = eta, beta = beta, theta = theta)
}
##################################################

