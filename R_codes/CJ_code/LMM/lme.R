##library(nlme)
##lme
##random effect model structure
if (k==1){
	remd <- pdIdent(~Z-1)
} else {
remd <- list()
m0 <- 0
for (i in 1:k){
	j <- (m0+1):(m0+d[i])
	remd[[i]] <- pdIdent(as.formula(paste0("~", paste0(colnames(Z)[j], collapse = "+"), "-1")))
	m0 <- m0 + d[i]
	}
remd <- pdBlocked(remd)
}

##one group variable
group <- rep("1", length(y))

fit <- lme(y ~ X-1, data = as.data.frame(G), random = list(group = remd), method = "REML")
	#summary(fit)
	#summary(fit)$sigma
	#summary(fit)$sigma^2
	theta.lme <- VarCorr(fit)[, "Variance"]
	theta.lme <- as.numeric(theta.lme[!duplicated(theta.lme)])

