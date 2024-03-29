##Use outputs from LMM fit to test fixed effects or contrasts of fixed effects.
##Arguments:
##- fit, outputs from lmmfit or lmmfitSS which contains 
##  'coef' (estimates of fixed effects), a matrix with rows representing the fixed effects and columns the different response variables in the model. 
##  'cov' (covariance matrix of the fixed effects), an array of three dimmesions for different response variables in the model. 
##  'df' (residual degree of freedom in the linear model).
##- index, integer or character vector indicating which fixed effects are to be tested. 
##  By default index consists of all of the fixed effects. Ignored if contrast is not NULL.
##- contrast, matrix with columns corresponding to contrasts of the fixed effects to be tested.
##- alternative, a character string specifying the alternative hypothesis, 
##  one of "two.sided", "greater" or "less".
##
lmmtest <- function(fit, index, contrast = NULL, alternative = c("two.sided", "less", "greater")){
alternative <- match.arg(alternative)
	if (is.null(contrast)){
		if (missing(index)) index <- 1:nrow(fit$coef)
		contrast <- diag(nrow(fit$coef))
		colnames(contrast) <- rownames(fit$coef)
		contrast <- contrast[, index, drop = FALSE]
	} 
	
	tval <- t(contrast)%*%fit$coef
	for (j in 1:ncol(fit$coef)){
		tval[, j] <- tval[, j]/sqrt(diag(t(contrast)%*%fit$cov[,,j]%*%contrast))
		}
		
	df <- fit$df
	if (alternative == "less") {
        pval <- pt(tval, df)
    } else if (alternative == "greater") {
        pval <- pt(tval, df, lower.tail = FALSE)
    } else pval <- 2 * pt(-abs(tval), df)

	rownames(tval) <- paste0(rownames(tval), "_t")
	rownames(pval) <- paste0(rownames(pval), "_pvalue")
      	
cbind(t(tval), t(pval))
}
