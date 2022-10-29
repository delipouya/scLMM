
##################################################
##Jacobi (diagonal) preconditioned conjugate gradients for a correlation matrix
##Q = A'A + lambda*I
##A'A = DG'RGD, a genotypic similarity matrix (similarity or correlation between genotypes)
##A'A = GDDG', a kinship matrix (genetic resemblance matrix between individuals)
##
##A = RGD for genotypic correlation (similarity)
##A = DG' for Kinship matrix
##
##Conjugate gradient method for computing inverse matrix multiplied by a vector:
## x = Q^{-1}b, or Qx = b, 
##x and b may be an n-by-p matrix.
##
##Jacobi (diagonal) preconditioning matrix 
## M = diag(Q) = diag(lambda + colSums A^2)
##################################################

pcgDiag <- function(A, b, x = 0*b, lambda = 1, M = lambda + colSums(A^2),
		i.rnew = 50, maxit = ncol(A), epsilon = 1e-8, verbose = FALSE){
##M = lambda + colSums(A^2): the diagonal of predconditioning matrix (Jacobi)
##M = 1, without precondioning.
	i <- 0
	r <- b - t(A)%*%(A%*%x) - lambda*x  #b-Qx, '%*%' turns to being as a matrix
	d <- sweep(r, MARGIN = 1, STATS = M, FUN = "/")  #M^{-1}r instead of d = r
	r20 <- colSums(r*d) #r^T*d instead of r^T*r
	r2new <- r20
	
	epsilon <- min(r20)*epsilon^2  #min(colSums(r^2))*epsilon^2
	while ((i < maxit) & (max(r2new) > epsilon)){
		Qd <- t(A)%*%(A%*%d) + lambda*d
		alpha <- r2new/(colSums(d*Qd))
		x <- x + sweep(d, MARGIN = 2, STATS = alpha, FUN = "*") #alpha*d
		if (i%%i.rnew == 0) r <- b - t(A)%*%(A%*%x) - lambda*x #b-Qx
		else r <- r - sweep(Qd, MARGIN = 2, STATS = alpha, FUN = "*")  #alpha*Qd
	
		s <- sweep(r, MARGIN = 1, STATS = M, FUN = "/") 
		r2old <- r2new
		r2new <- colSums(r*s) #instead of r^2
		beta <- r2new/r2old
		d <- s + sweep(d, MARGIN = 2, STATS = beta, FUN = "*")  #instead of r + beta*d
		
		i <- i + 1
	}

if (verbose) cat(paste0("rsquare*epsilon = ", epsilon, ". Number of iterations = ", i, "."))
return(x)
}

##################################################
##source("R/pcgDiag.R")
if (FALSE) {
	A <- matrix(1:16, 4,4)
	A <- matrix(rnorm(4*6, m=0, s=5), 4, 6)
	b <- 1:ncol(A)

	Q <- t(A)%*%A
	M <- Q%*%Q - Q
	D <- diag(1/(1 + rowSums(A^2)))
	M <- - t(A)%*%D%*%A
	
	lambda <- 50
	lambda <- 0.1
	diag(Q) <- diag(Q) + lambda
	diag(M) <- diag(M) + lambda
	
	#svd(Q)
	eigen(Q)
	eigen(M%*%Q)
	
x0 <- solve(Q)%*%b
x1 <- pcgDiag(A, b, lambda = lambda, verbose = T)
x2 <- pcgDiag(A, b, lambda = lambda, M = 1, verbose = T)
	abs(x1-x0)
	abs(x2-x0)
}
##################################################