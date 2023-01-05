# Load the necessary libraries
library(lme4)
library(ggplot2)

# Generate some fake data
set.seed(123)
n <- 10 # number of students
x <- rnorm(n*10) # predictor variable
student <- (rep(1:n, each=10)) # student ID
y <- x + rnorm(n*10, sd=0.5) + 0.5*student # response variable

student<- factor(student)
# Fit the mixed model
fit <- lmer(y ~ x + (1|student), data=data.frame(x, y, student))

# Extract the fixed and random effects
fe <- fixef(fit)
re <- ranef(fit)

# Calculate predicted values for each student
pred <- matrix(NA, nrow=n, ncol=2)
colnames(pred) <- c("intercept", "slope")
for (i in 1:n) {
  pred[i,1] <- fe[1] + re$student[i,]
  pred[i,2] <- fe[2]
}

df = data.frame(x, y, student)
# Plot the results
ggplot(df, aes(x=x, y=y, color=student)) +
  geom_point() +
  geom_abline(intercept=fe[1], slope=fe[2], linetype="dashed") +
  geom_abline(aes(intercept=intercept, slope=slope), data=as.data.frame(pred), linetype="dotted")+
  ggtitle('student ID as random effect inn LMM')

  


set.seed(123)
n <- 10 # number of students
x <- rnorm(n*10) # predictor variable
student <- (rep(1:n, each=10)) # student ID
y <- x + rnorm(n*10, sd=0.5) + 0.5*student # response variable

student<- factor(student)

# Fit the linear model
fit <- lm(y ~ x + student, data=data.frame(x, y, student))

# Extract the coefficients
coef <- coef(fit)

# Calculate predicted values for each student
pred <- matrix(NA, nrow=n, ncol=2)
for (i in 1:n) {
  pred[i,1] <- coef[1] + coef[i+2]
  pred[i,2] <- coef[2] + coef[i+2]
}

# Plot the results
ggplot(data.frame(x, y, student), aes(x=x, y=y, color=student)) +
  geom_point() +
  geom_abline(intercept=coef[1], slope=coef[2], linetype="dashed") +
  geom_abline(aes(intercept=intercept, slope=slope), data=pred, linetype="dotted")

