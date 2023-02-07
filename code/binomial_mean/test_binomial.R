library(vebpm)
library(fastGHQuad)
library(ebnm)
library(Rfast)
library(smashrgen)
set.seed(12345)
n=1024
mu = c(rep(-1,n/2),rep(1,n/4),rep(-1,n/4))
p = 1/(1+exp(-mu))
nb = rep(100,n)
x = rbinom(n,nb,p)
plot(x/nb,col='grey80')
lines(p,col='grey60')
fit = binom_smooth_split(x,nb,m_init = logit((x+1)/(nb+2)),sigma2_init=1,tol=1e-8,maxiter = 100)
#fit = binom_smooth_split(x,nb,m_init = logit(x/nb),sigma2_init=0.5,tol=1e-8)
#fit = binomial_mean_splitting(x,nb,printevery = 1)
lines(fit$posterior$mean_smooth)


xx = (1:n)/n
d <- data.frame(xx,x) ## not absolutely necessary but good practice
library(mgcv)
m1 <- gam(x~s(xx),family="binomial",data=d)
lines(m1$fitted.values)
