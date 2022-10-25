
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

n = 1000
mu = c(rep(0,n*0.8),rnorm(n*0.2,0,3))
#p = exp(mu)/(1+exp(mu))
#r = 10
#x = rnbinom(n,r,1-p)
plot(mu,col='grey80')
plot(x,col='grey80')
#x = rpois(n,exp(mu))


n = 1000
mu = c(rep(2,n*0.8),rnorm(n*0.2,2,1))
#p = exp(mu)/(1+exp(mu))
#r = 10
#x = rnbinom(n,r,1-p)
x = rpois(n,exp(mu))
plot(mu,col='grey80')
plot(x,col='grey80')
#

temp = pois_mean_penalized_inversion(x)
plot(mu,col='grey80')
lines(temp$posteriorMean)

temp = pois_mean_penalized_compound(x)

temp = ash_pois(x,link='log')
plot(mu,col='grey80')
lines(log(temp$result$PosteriorMean))

temp = ash_pois(x,link='identity')
plot(mu,col='grey80')
lines(log(temp$result$PosteriorMean))

temp = pois_mean_split(x)
temp$ebnm_res$fitted_g
plot(mu,col='grey80')
lines(temp$posteriorMean)

temp = pois_mean_GMGM(x)
plot(mu,col='grey80')
lines(temp$posteriorMean)

temp = pois_mean_GG(x)
plot(temp$obj_value,type='l')
plot(mu,col='grey80')
lines(temp$posteriorMean)

temp = pois_mean_GMG(x)
plot(temp$obj_value,type='l')
plot(mu,col='grey80')
lines(temp$posteriorMean)

temp = pois_mean_log1exp(x)
plot(temp$obj_value,type='l')
plot(mu,col='grey80')
lines(log(log(1+exp(temp$posteriorMean))))

temp = nb_mean_polya_gamma(x,r=100,ebnm_params = list(prior_family = "point_laplace"))
plot(temp$obj_value,type='l')
temp$ebnm_res$fitted_g
plot(mu,col='grey80')
lines(temp$poisson_log_mean_est)

temp = nb_mean_lower_bound(x,r=100,ebnm_params = list(prior_family = "point_laplace"))
plot(temp$obj_value,type='l')
temp$ebnm_res$fitted_g
plot(mu,col='grey80')
lines(temp$poisson_log_mean_est)

temp = nb_ash_pg(x,r=100)
plot(temp$ELBO,type='l')
plot(mu,col='grey80')
#plot(exp(mu),col='grey80')
#lines(temp$poisson_mean_est)
lines(temp$poisson_log_mean_est)


plot(mu,col='grey80',ylim=range(c(mu,temp$poisson_log_mean_est)))
lines(temp$poisson_log_mean_est)
plot(temp$obj_value,type='l')

