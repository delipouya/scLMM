
######## Evaluating the residuals ########
list_of_residuals = lapply(models, function(x)(residuals(x)))
residual_df = t(do.call(cbind, list_of_residuals))
residual_df[1:10, 1:10]
hist(residual_df, 100)

