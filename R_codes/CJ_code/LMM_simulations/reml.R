##################################################
##Iterative REML algorithms: 
##1) Newton-Raphson (NR)
##2) Fisher's information score (FS)
##3) Average information (AI)
##################################################

reml <- function(y, X, Z, d, s = rep(0.5, length(d) + 1), 
	method = c("FS", "AI", "NR"), niter = 5)
{
##s = (s1, ...,sk, s_{k+1}), 
##si = sigma_i^2, s_{k+1} = sigma_0^2
##d = (m1,...,mk), mi = ncol(Zi), m1+...+mk = ncol(Z)
##method: 1) Newton–Raphson (NR), 2) Fisher scoring (FS), 3) Average information (AI)

method <- match.arg(method)	
k <- length(d)  #number of regions
n <- length(y)  
p <- ncol(X)

	
##########
##estimated variance components in iterations
theta.iter <- matrix(NA, nrow = k+1, ncol = niter)
rownames(theta.iter) <- paste0("theta", 1:(k+1))
colnames(theta.iter) <- 1:niter

##first derivatives of log likelihood in iterations
dl.iter <- matrix(NA, nrow = k+1, ncol = niter)
	
for (iter in 1:niter){
			
	V <- diag(rep(s[k+1], n))
	mi <- 0
	for (i in 1:k){	
		Zi <- Z[, (mi+1):(mi+d[i])]
		V <- V + s[i]*Zi%*%t(Zi)
		mi <- mi + d[i]
		}
	Vin <- solve(V)
	Rs <- Vin - (Vin%*%X)%*%solve(t(X)%*%Vin%*%X)%*%(t(X)%*%Vin)
	Ry <- Rs%*%y
		
	##########
	##computing matrix H, I, and Ia
	##S = RθViRθVj
	##Sy = y'RθViRθVjRθy 
	##I = S/2, Ia = Sy/2, H = S/2 - Sy 
	S <- matrix(NA, k+1, k+1)	
	Sy <- matrix(NA, k+1, k+1)	
	dl <- rep(NA, k+1)
	
	mi <- 0
	for (i in 1:k){	
	Zi <- Z[, (mi+1):(mi+d[i])]
	Vi <- Zi%*%t(Zi)
	RVi = Rs%*%Vi
		
	##dl
	dl[i] <- (t(Ry)%*%Vi%*%Ry - sum(diag(RVi)))/2
	
	mj <- 0
	for (j in 1:i){
		Zj <- Z[, (mj+1):(mj+d[j])]
		Vj <- Zj%*%t(Zj)
		RVj <- Rs%*%Vj
		
		S[i, j] <- sum(diag(RVi%*%RVj))/2
		S[j, i] <- S[i, j]
		##Sy
		Sy[i, j] <- t(Ry)%*%Vi%*%Rs%*%Vj%*%Ry/2
		Sy[j, i] <- Sy[i, j]
		
		mj <- mj + d[j]
		}
		
	j <- k+1
		S[i, j] <- sum(diag(RVi%*%Rs))/2
		S[j, i] <- S[i, j]
		##Sy
		Sy[i, j] <- t(Ry)%*%Vi%*%Rs%*%Ry/2
		Sy[j, i] <- Sy[i, j]

	mi <- mi + d[i]
	}
	
	i <- k+1
	S[i, i] <- sum(diag(Rs%*%Rs))/2
	##Sy
	Sy[i, i] <- t(Ry)%*%Rs%*%Ry/2
	##dl
	dl[i] <- (t(Ry)%*%Ry - sum(diag(Rs)))/2
	
	##########
	H = S - 2*Sy
	FS = S
	AI = Sy
	
	if (method == "NR") {
		M <- -H
	} else {
	if (method == "AI") M <- AI
	else M <- FS
	}
	
	##The H and AI matrices can be singular.
	#Minv <- solve(M)
	Minv <- ginv(M)
	
	s <- s + Minv%*%dl
	theta.iter[, iter] <- s
	dl.iter[, iter] <- dl
	}

list(dl.iter = dl.iter, theta.iter = theta.iter, H = H, FS = FS, AI = AI, theta = theta.iter[, niter], method = method)
}

##################################################
##source("R/reml.R")
##out <- reml(y, X, Z, d)
##out <- reml(y, X, Z = G, d, niter = 10)

##################################################
##notes

if (FALSE) {

"https://www.animalgenome.org/bioinfo/resources/manuals/ASReml2/"
"https://www.animalgenome.org/bioinfo/resources/manuals/ASReml2/htmlhelp/asreml/critical.htm"

"Singularities in Average Information matrix
In release 1.00, if singularities were present in the AI matrix, a generalised inverse was used which effectively conditioned on whichever parameters were identified as singular. 
ASReml now aborts processing if singularities appear unless the !AISINGULARITIES qualifier is set. 
Which particular parameter is singular is reported in the variance component table printed in the .asr file."

"The most common reason for singularities is that the user has overspecified the model and is likely to misinterpret the results if not fully aware of the situation. Overspecification will occur in a direct product of two unconstrained variance matrices, when a random term is confounded with a fixed term and when there is no information in the data on a particular component. The best action is to reform the variance model so that the ambiguity is removed, or to fix one of the parameters in the variance model so that the model can be fitted. Only rarely will it be reasonable to specifiy the !AISINGULARITIES qualifier."

"Mary J. Lindstrom and Douglas M. Bates (1988): Unlike the EM algorithm, the NR algorithm is not guaranteed to converge. The NR algo- rithm, however, implemented as we suggest (optimizing the profile likelihood with respect to , and the entries of the Cholesky factor of D), will converge to a local maxi- mum of the likelihood surface very consistently."

"Liu, Chapter 4: Sometimes, however, the researcher might encounter failure of convergence in the iterative processes or the occurrence of convergence at unrealistic parameter val- ues when performing the NR algorithm. "

} ##FALSE

