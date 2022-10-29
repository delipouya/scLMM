
##################################################
library(lme4)
library(nlme)
library(MASS)


##R functions
source("simulate.covar.geno.R")
source("reml.R")

##################################################
##Simulate data: y, X, G
##  X: covariates
##  beta: fixed effects
##  G: genotypes, G = [G1, ... , G_p], Gi is a matrix of genotypes in a region (kernel).
##  y: phenotypes, y = X*beta + G*b + e
##

SEED <- round(runif(1, 1e4, 1e5))
set.seed(SEED)

#####
##simulate X, beta, and G using function simulate.covar.geno
##inputs for simulate.covar.geno
n <- 1000       #sample size
p <- 1         #number of covariates
d <- c(2, 3)  #numbers of SNPs in each kernels #len(d)=number of kernels
rho <- 0.2#0.5     #correlation between SNPs in each kernel (region)

XG <- simulate.covar.geno(n = n, p = p, d = d, rho = rho, fixed.m = 0, fixed.sd = 0)
X <- XG$X
beta <- XG$beta
G <- XG$G

head(G)
summary(data.frame(G))
head(X)
#sum(G == 0)/length(G)
#sum(G == 1)/length(G)
#sum(G == 2)/length(G)


#####
##simulate y
sigma <- c(0.4, 0.6, 1)  #square root of variance components
sigma^2  #true variance components

##random effects
b <- sapply(1:length(d), function(i) rnorm(d[i], sd = sigma[i]))
b <- unlist(b)  
##phenotypes
y_mat = matrix(data = NA, nrow = 200, ncol = n)
y <- X%*%beta + G%*%b + rnorm(n, sd = sigma[length(d) + 1])

write.csv(y, 'y.csv')
write.csv(X, 'X.csv')
write.csv(G, 'G.csv')


rm(XG, b)
nrow

##################################################
##compare three reml iterative algorithms with lme{nlme}
##Iterative algorithms: 
##1) Newton-Raphson (NR)
##2) Fisher's information score (FS)
##3) Average information (AI)
	
	
#####
##use lme function in nlme package

##random effect model structure
k <- length(d) ##number of kernels (regions)

if (k==1){
	lme.md <- pdIdent(~G-1)
} else {
	lme.md <- list()
	m0 <- 0
	for (i in 1:k){
	j <- (m0+1):(m0+d[i])
	lme.md[[i]] <- pdIdent(as.formula(paste0("~", paste0(colnames(G)[j], collapse = "+"), "-1")))
	m0 <- m0 + d[i]
	}
lme.md <- pdBlocked(lme.md)
}

##lme needs that sample size is greater than the number of random effect variables.
##n >= sum(d)

group <- rep("1", length(y)) #one group variable
fit <- lme(y ~ X-1, data = as.data.frame(G), random = list(group = lme.md), method = "REML")
	#summary(fit)
	#summary(fit)$sigma
	#summary(fit)$sigma^2
	v <- VarCorr(fit)[, "Variance"]
	theta.lme <- as.numeric(v[!duplicated(v)])


#####
##REML
##use R function reml

k <- length(d)# number of kernels (regions)
SEED.s0 <- round(runif(1, 1, 1e4))
set.seed(SEED.s0)
s0 <- runif(k + 1)  #initial values of the variance components
niter <- 15 #20 #number of iterations


##MINQE
#minqe.method <- "IU" 
#mq <- minqeRao(y, X, G, d = d , s0 = s0, method = minqe.method, niter = niter)
	
##REML
par(mfrow = c(3,2), mar = c(4.5, 4.5, 2.5, 2))

for (iter.method in c("FS", "AI", "NR")){
	ml <- reml(y, X, Z = G, d, s = s0, method = iter.method, niter = niter)
	
	u <- ml$theta.iter
	ylim <- range(sigma^2, u)
	plot(u[1, ], ylim = ylim, type = "n", xlab = "Number of iterations", ylab = "Variance components", main = paste0("REML: ", iter.method, " (true in grey, lme in blue)"))
	for (i in 1:nrow(u)) {
		lines(u[i, ], lty = 1)
		points(u[i, ], pch =16, cex = 0.7)
		}
	points(rep(1,nrow(u)), sigma^2, pch = 8, col = "grey")
	points(rep(1,nrow(u)), theta.lme, pch = 1, col = "blue")
		
	u <- ml$dl
	plot(u[1, ], ylim = range(u), type = "n", xlab = "Number of iterations", ylab = "First derivatives", main = paste0("REML: ", iter.method))
	for (i in 1:nrow(u)) {
		lines(u[i, ], lty = i)
		points(u[i, ], pch =16, cex = 0.7)
		}
	}
	

	
##################################################
##Multiple simulations for REML with FS method
##- Design matrices X and G, fixed effects beta, and variance components are the same for all simulations.
##- Phenotypes y are varied based on the random effects and residual errors.


sigma <- c(0.4, 0.6, 1)  #square root of variance components
sigma^2 #true variance components

##simulate X, beta, and G using function simulate.covar.geno
##inputs for simulate.covar.geno
n <- 100       #sample size
p <- 1         #number of covariates
d <- c(5, 10)  #numbers of SNPs in each kernels 
rho <- 0.5     #correlation between SNPs in each kernel (region)

XG <- simulate.covar.geno(n = n, p = p, d = d, rho = rho, fixed.m = 0, fixed.sd = 0)
	X <- XG$X
	beta <- XG$beta
	G <- XG$G
	
	
##use lme function in nlme package
##random effect model structure
k <- length(d) ##number of kernels (regions)
if (k==1){
	lme.md <- pdIdent(~G-1)
} else {
	lme.md <- list()
	m0 <- 0
	for (i in 1:k){
	j <- (m0+1):(m0+d[i])
	lme.md[[i]] <- pdIdent(as.formula(paste0("~", paste0(colnames(G)[j], collapse = "+"), "-1")))
	m0 <- m0 + d[i]
	}
lme.md <- pdBlocked(lme.md)
}
	
	
##	
nsimu <- 100 #number of simulation
iter.method <- "FS"
niter <- 8

theta <- array(dim = c(length(d) + 1, niter, nsimu)) 
dl <- array(dim = c(length(d) + 1, niter, nsimu)) 
theta.lme <- matrix(nrow = length(d) + 1, ncol = nsimu)

#####
for (nrep in 1:nsimu) {
	##simulate y
	##random effects
	b <- sapply(1:length(d), function(i) rnorm(d[i], sd = sigma[i]))
	b <- unlist(b)  
	##phenotypes
	y <- X%*%beta + G%*%b + rnorm(n, sd = sigma[length(d) + 1])
	rm(b)
	
	##REML
	ml <- reml(y, X, Z = G, d, s = s0, method = iter.method, niter = niter)
	theta[,,nrep] <- ml$theta.iter
	dl[,,nrep] <- ml$dl
	
	##lme
	group <- rep("1", length(y)) #one group variable
	fit <- lme(y ~ X-1, data = as.data.frame(G), random = list(group = lme.md), method = "REML")
	v <- VarCorr(fit)[, "Variance"]
	theta.lme[, nrep] <- as.numeric(v[!duplicated(v)])
}

##
u <- apply(theta, 1:2, mean)
#v <- apply(theta, 1:2, sd)
v <- apply(theta, 1:2, sd)/sqrt(nsimu)
u.lme <- apply(theta.lme, 1, mean)

par(mfrow = c(nrow(u), 1), mar = c(5,4,2,2))
for (i in 1:nrow(u)) {
	ylim <- range(u[i, ], u[i, ] + 1.96*v[i, ], u[i, ] - 1.96*v[i, ])
	plot(u[i, ], ylim = ylim, xlab = "Number of iterations", ylab = paste0("Variance ", i), type = "l")
	points(u[i, ], pch = 16)
	lines(u[i, ] + 1.96*v[i, ], lty = 2)
	lines(u[i, ] - 1.96*v[i, ], lty = 2)	

	points(1, sigma[i]^2, pch = 8, col = "red", cex = 1.5)
	points(1, u.lme[i], pch = 1, col = "blue", cex = 1.5)
	
	mtext("True value in red, lme in blue. Dashed lines represent 95% CI of the REML estimates.", cex = 0.9)
	}
		

##################################################