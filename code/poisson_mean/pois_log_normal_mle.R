library(fastGHQuad)

pois_log_normal_mle = function(y,prior_mean,prior_var_init = 1,n_gh=20){
  gh_points = gaussHermiteData(n_gh)
  res = optim(log(prior_var_init),
              pois_log_normal_mle_obj,
              y=y,
              prior_mean=prior_mean,
              gh_points = gh_points,
              method = 'L-BFGS-B')
  exp(res$par)
}


pois_log_normal_mle_obj = function(prior_var_log,y,prior_mean,gh_points){
  prior_var = exp(prior_var_log)
  K = length(gh_points$x)
  n = length(y)
  const = -sum(lfactorial(y)) - n*log(pi)/2
  temp = (sqrt(2*prior_var)*tcrossprod(y,gh_points$x) + matrix(y*prior_mean,nrow=n,ncol=K,byrow=F)-exp(sqrt(2*prior_var)*matrix(gh_points$x,nrow=n,ncol=K,byrow=T) + matrix(prior_mean,nrow=n,ncol=K,byrow=F)))
  obj = const + sum(log(exp(temp)%*%gh_points$w))
  return(-obj)
}

set.seed(12345)
n = 10000
#mu = rnorm(n,-6,1)
mu = l0 + f0[1]
n = length(mu)
sigma2 = 0
gh_points = gaussHermiteData(50)
y = rpois(n,exp(rnorm(n,mu,sqrt(sigma2))))
sum(y==1)
res2 = vebpm::ebpm_normal(y,g_init=list(mean=mu,var=NULL),fix_g=c(T,F))
res2$fitted_g$var

pois_log_normal_mle(y,mu)


objs = c()
prior_vars = seq(1e-8,4,length.out=1000)
for(i in 1:1000){
  objs[i] = pois_log_normal_mle_obj(log(prior_vars[i]),y,mu,gh_points)
}
plot(prior_vars,objs,type = 'l')
