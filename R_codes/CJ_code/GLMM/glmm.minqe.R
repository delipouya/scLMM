
##################################################
##MINQE computing from simplified Rao's formulas with 
##- (1) trace estimation using R function minqeGeno.tr
##- (2.1) conjugate gradients using R function cg.mat 
##- (2.2) Jacobi (diagonal) preconditioned conjugate gradients using R function pcgDiag 

glmm.minqe <- function(y, X, G, d, family = c("binomial", "poisson"), phi = NULL, method = c("IU", "I"), 
		niter = 10, niter.minqe = 5, nsamples = 100, sampleDist = rnorm, epsilon = 1e-8, plot.it = FALSE) {
##d = (m1,...,mk), numbers of geneotypes
##Using Rademacher distribution (Hutchinson's estimator) has lower variance compared to the Gaussian.
##Trace estimation converges fast as total number of genotypes becomes larger.
##Hutchinson's estimator:
##Rademacher distribution = 2*Bernoulli - 1
##sampleDist <- function(n) 2*rbinom(n, size = 1, prob = 0.5) - 1
##minqeGeno
##

family <- match.arg(family)

##########
##functions:
## inverse link
## weights = {1/[a_iv(mu_i)g'(mu_i)^2]}
## working vector = g(mu) + g'(mu)(y-mu)
##########

if (family == "binomial"){
	size <- y$size
	y <- y$y
	ginv <- function(eta) size/(1 + exp(-eta))
	weight <- function(mu) mu*(1-mu/size)
	working <- function(mu) log(mu/(size - mu)) + size*(y - mu)/(mu*(size - mu))
	mu <- (y+0.5*size)/2
}

if (family == "poisson"){
	ginv <- function(eta) exp(eta)
	weight <- function(mu) mu
	working <- function(mu) log(mu) + (y - mu)/mu
	mu <- y+0.5
}

##initial values
theta = rep(1, length(d)+1)

##########
beta <- matrix(NA, nrow = niter, ncol = ncol(X))
thetait <- matrix(NA, nrow = niter, ncol = length(theta)) 
for (i in 1:niter){
	w <- sqrt(weight(mu))
	zw <- w*working(mu)
	Xw <- sweep(X, MARGIN = 1, STATS = w, FUN = "*")
	Gw <- sweep(G, MARGIN = 1, STATS = w, FUN = "*")
	
	lmmfit <- minqeGeno(zw, Xw, Gw, d = d, theta = theta, phi = phi, method = method, niter = niter.minqe,
			nsamples = nsamples, sampleDist = sampleDist, epsilon = epsilon)
	beta[i, ] <- lmmfit$beta
	theta <- lmmfit$theta
	thetait[i, ] <- theta
	eta <- sweep(lmmfit$eta, MARGIN = 1, STATS = w, FUN = "/")
	mu <- ginv(eta)
	}


if (plot.it){
par(mfrow = c(2, 1), mar = c(4.1, 4.1, 3.1, 2.1))
plot(beta[,1], ylim = range(beta), xlab = "Iterations", ylab = "beta", type = "n")
cols <- rainbow(max(7, ncol(beta)))
for (j in 1:ncol(beta)) {
	lines(beta[, j], col = cols[j])
	points(beta[, j], col = cols[j], pch = 16, cex = 0.8)
	}

plot(thetait[,1], ylim = range(thetait), xlab = "Iterations", ylab = "theta", type = "n")
cols <- rainbow(max(7, ncol(thetait)))
for (j in 1:ncol(thetait)) {
	lines(thetait[, j], col = cols[j])
	points(thetait[, j], col = cols[j], pch = 16, cex=0.8)
	}
}
	
	
list(beta = beta[niter, ], theta = theta)
}
##################################################

