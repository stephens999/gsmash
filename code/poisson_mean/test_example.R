set.seed(12345)
n = 1000
w = 0.2
mu = log(c(rgamma(n*(1-w),2,2),rgamma(n*w,10,2)))
x = rpois(n,exp(mu))

set.seed(12345)
lambda = c(rep(0.1,n*(1-w)),rep(10,n*w))
x = rpois(n,lambda)
out = nb_mean_polya_gamma(x,r=1000,maxiter = 10000)
plot(out$obj,type='l')
plot(x,main='observation',col='grey80',pch=20)
lines(lambda,col='grey50')
lines(out$mean_est,col=4)



out = pois_mean_log1exp(x,maxiter = 300)
out = pois_mean_GMGM(x,tol=1e-5,maxiter = 500)


#out = pois_mean_GG(x,tol=1e-5)
plot(out$obj,type='l')

out$beta
out$sigma2

par(mfrow=c(3,1))
plot(x,main='observation',col='grey80',pch=20)
plot(mu,main='parameter',col='grey80',pach=20)
plot(out$m,main='posterior mean',col='grey80',pch=20)

plot(log(x))
lines(out$m)
lines(mu,col=3)


out2 = pois_mean_GMG(x,tol=1e-3)


