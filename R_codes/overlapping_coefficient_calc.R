min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

#CASE 1
mu1 <- 0
mu2 <- 0
sd1 <- .1
sd2 <- .15

xs <- seq(min(mu1 - 3*sd1, mu2 - 3*sd2), max(mu1 + 3*sd1, mu2 + 3*sd2), .01)
f1 <- dnorm(xs, mean=mu1, sd=sd1)
f2 <- dnorm(xs, mean=mu2, sd=sd2)
samples1 = xs

plot(xs, f1, type='l', ylim=c(0, max(f1,f2)), ylab='density', xlab='Factor score',font=10 )
lines(xs, f2, lty='dotted')
ys <- min.f1f2(xs, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])
polygon(xs, ys, col='gray')

n = 20000
f1_sample <- rnorm(n, mean=mu1, sd=sd1)
f2_sample <- rnorm(n, mean=mu2, sd=sd2)
df = data.frame(sample=c(f1_sample,f2_sample),
                label=c(rep('Covariate-1', n), rep('Covariate-2', n)))
ggplot(df, aes(x=sample,fill=label))+geom_density(alpha=0.5)+theme_classic()+
  scale_fill_manual(values = c('orange', 'maroon'))+theme(text = element_text(size=20))+xlab('Factor score')





print(paste('OVL:',integrate(min.f1f2, -Inf, Inf, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)$value,sep=))
print(paste('FN:',abs(sd1^2 - sd2^2),sep=''))

#CASE 2
mu1 <- 0
mu2 <- 0
sd1 <- .1
sd2 <- .05
sample2 = xs
## new case
mu1 <- 0.45
mu2 <- 1
sd1 <- .25
sd2 <- .1




n = 30000
f1_sample <- rnorm(n, mean=mu1, sd=sd1)
f2_sample <- rnorm(n, mean=mu2, sd=sd2)
df = data.frame(sample=c(f1_sample,f2_sample),
                label=c(rep('Covariate-1', n), rep('Covariate-2', n)))
ggplot(df, aes(x=sample,fill=label))+geom_density(alpha=0.5)+theme_classic()+
  scale_fill_manual(values = c('orange', 'maroon'))+theme(text = element_text(size=20))+xlab('Factor score')



plot(sample2, samples1)
xs <- seq(min(mu1 - 3*sd1, mu2 - 3*sd2), max(mu1 + 3*sd1, mu2 + 3*sd2), .01)
f1 <- dnorm(xs, mean=mu1, sd=sd1)
f2 <- dnorm(xs, mean=mu2, sd=sd2)

plot(xs, f1, type='l', ylim=c(0, max(f1,f2)), ylab='density', xlab="Factor score")
lines(xs, f2, lty='dotted')
ys <- min.f1f2(xs, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])
polygon(xs, ys, col='gray')
sample3 = xs

print(paste(“OVL:',integrate(min.f1f2, -Inf, Inf, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)$value,sep=' “))
print(paste(“FN:',abs(sd1^2 – sd2^2),sep=' “))

#Compare over range of SD values
x = matrix(seq(from=.01,to=.2,by=.01))
fn = apply(x,1,function(s1) abs(.01 – s1^2))
ovl = apply(x,1,function(s1) integrate(min.f1f2, -Inf, Inf, mu1=0, mu2=0, sd1=s1, sd2=0.1)$value)

toGraph = cbind(x,fn,ovl)

#indiffernce compared to SD=0.14
fn = abs(0.01 – .14^2)
ovl = integrate(min.f1f2, -Inf, Inf, mu1=0, mu2=0, sd1=.14, sd2=0.1)$value

plot(toGraph[,2],toGraph[,3],xlab='FN',ylab='OVL',ylim=c(0,1.1))
abline(h=ovl,lty='dotted')
abline(v=fn,lty='dotted')
text(toGraph[,2],toGraph[,3],toGraph[,1],col='red',pos=1)