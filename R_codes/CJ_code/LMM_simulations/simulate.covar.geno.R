##################################################
##Generate data
##X and beta: covariates and cofficients (fixed effects)
##G: genotypes
##################################################

simulate.covar.geno <- function(n = 25, p = 1, d = 5, fixed.m = 1, fixed.sd = 0.1, rho = 0.5, p0 = 0.7, p1 = 0.2){
##d = (q1,...,qk), q=q1+...+qk
##p0+p1+p2=1
##p0=Pr(SNP=0)
##p1=Pr(SNP=1)
##p2=Pr(SNP=2)
##rho = cor(SNP_i, SNP_{i+1})
##sigma.squared: variance components
##X = [X1,...,X_p]

k <- length(d)

##design matrix X and fixed effects beta
X <- matrix(rnorm(n*p), n, p)
colnames(X) <- paste0("X", 1:p)
beta <- rnorm(p, m = fixed.m, sd = fixed.sd)

##design matrix G and random Effects b, 
##correlated G within groups
##rho <- 0.5

G <- NULL
for (j in 1:k){
	A <- diag(d[j])
	for (i in 1:min(9, (d[j]-1))){
	A[row(A) == col(A) - i] <- rho^i
	A[row(A) == col(A) + i] <- rho^i
	}
	B <- chol(A)
	#range(t(B)%*%B-A)

	G0 <- t(t(B)%*%matrix(rnorm(n*d[j]), d[j], n))
	#cor(G0)
	qnt <- quantile(G0, probs = c(p0, p0+p1))
	Gj <- matrix(2, n, d[j])
	Gj[G0 <= qnt[2]] <- 1
	Gj[G0 <= qnt[1]] <- 0	
	#cor(Gj)
	rm(G0)
	
	G <- cbind(G, Gj)
	rm(Gj)
	}
colnames(G) <- paste0("SNP", 1:ncol(G))

list(X = X, beta = beta, G = G, d = d)
}

##################################################
##geno = simulate.covar.geno()
##geno = simulate.covar.geno(n=10, p=2, d=c(2, 5))
##	sum(geno$G>0)/(length(geno$G))
##	sum(geno$G==1)/(length(geno$G))
##	sum(geno$G==2)/(length(geno$G))
##