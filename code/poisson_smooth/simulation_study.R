#'@title generate Poisson sequence data
#'@param n_simu number of simulation reps
#'@param n length of Poisson seq
#'@param snr signal  to noise ratio
#'@param count_size max of exp(b)
#'@param smooth_func blocks, bumps, heavi,doppler
#'@return a list of
#'  \item{X:}{data}
#'  \item{L:}{mean}
#'  \item{other imputs:}{sigma2,snr, ...}
#'@import wavethresh
#'@import smashrgen
#'@import vebpm
#'@export
sim_data_smooth = function(n_simu,n=2^9,snr=3,count_size,smooth_func='blocks',seed=12345){
  set.seed(seed)
  b = DJ.EX(n=n,signal=1,noisy=FALSE,plotfn = FALSE)[[smooth_func]]
  b = b - min(b)
  b = b/(max(b)/log(count_size))
  sigma2 = var(b)/snr
  X = matrix(nrow=n_simu,ncol=n)
  L = matrix(nrow=n_simu,ncol=n)
  for(i in 1:n_simu){
    l = exp(b+rnorm(n,0,sd=sqrt(sigma2)))
    L[i,] = l
    X[i,] = rpois(n,l)
  }
  return(list(X=X,L=L,b=b,snr=snr,sigma2=sigma2,count_size=count_size,smooth_func=smooth_func,seed=seed))
}

#'@title compare methods for smoothing Poisson sequence
#'@export
simu_study_poisson_smooth = function(simdata,save_data=TRUE,
                                     method_list=c('vst','lik_exp','split_dwt',
                                                   'split_ndwt','smash'),
                                     n_cores = 10,
                                     filter.number = 1,
                                     family='DaubExPhase'){
  n_simu = nrow(simdata$X)
  n = ncol(simdata$X)
  n_method = length(method_list)
  res = mclapply(1:n_simu,function(i){
    fitted_model = vector('list',n_method)
    names(fitted_model) = method_list

    if('vst'%in%method_list){
      res_vst = try(smash_gen_pois(simdata$X[i,],transformation='vst',method='smash',
                                   filter.number = filter.number,family = family,maxiter=1))
      fitted_model$vst = res_vst
    }
    if('lik_exp'%in%method_list){
      res_lik = try(smash_gen_pois(simdata$X[i,],transformation='lik_expnasion',method='smash',
                                   filter.number = filter.number,family = family,maxiter=1))
      fitted_model$lik_exp = res_lik
    }
    if('split_dwt'%in%method_list){
      res_split = try(pois_smooth_split(simdata$X[i,],wave_trans='dwt',
                                   filter.number = filter.number,family = family))
      fitted_model$split_dwt = res_split
    }
    if('split_ndwt'%in%method_list){
      res_ndwt = try(pois_smooth_split(simdata$X[i,],wave_trans='ndwt',
                                        filter.number = filter.number,family = family))
      fitted_model$split_ndwt = res_ndwt
    }
    mse_smooth = simplify2array(lapply(fitted_model, function(x) {
      mse(x$posterior$mean_smooth, exp(simdata$b))}))
    mae_smooth = simplify2array(lapply(fitted_model, function(x) {
      mae(x$posterior$mean_smooth, exp(simdata$b))}))
    mse_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      mse(x$posterior$mean_latent_smooth,simdata$b)}))
    mae_latent_smooth = simplify2array(lapply(fitted_model, function(x) {
      mae(x$posterior$mean_latent_smooth,simdata$b)}))
    return(list(fitted_model=fitted_model,mse_smooth=mse_smooth,mae_smooth=mae_smooth,
                mse_latent_smooth=mse_latent_smooth,mae_latent_smooth=mae_latent_smooth))
  },mc.cores = n_cores
  )
  if(save_data){
    return(list(sim_data = sim_data, output = res))
  }else{
    return(res)
  }
}

#
# library(wavethresh)
# n = 512
# signal = 1
# y_list = DJ.EX(n=n,signal=signal,noisy=FALSE,plotfn = FALSE)
# b = y_list$blocks
# b = b - min(b) + 1
# set.seed(12345)
# snr = 1
# plot(b,type='l')
# sigma2 = var(b)/snr
# mu = b + rnorm(n,0,sqrt(sigma2))
#
# plot(mu,type='l')
# plot(exp(mu),type='l')
# lambda = exp(mu)
# y = rpois(n,lambda)
#
# library(gamlss)
# datax = data.frame(y=y,x=seq(0,1,length.out = n))
# mod = quote(gamlss(y~lo(~x,df=p[1]),sigma.formula = ~lo(~x,df=p[2]),data=datax,family = PO))
# op=find.hyper(model = mod, parameters  = c(1,1),lower=c(1,1),steps=c(0.1,0.1))
#
# fit = gamlss(y~lo(~x, df=op$par[1]),sigma.formula = ~lo(~x,df=op$par[2]),data=datax,family = PO())
# plot(y,col='grey80')
# lines(exp(b),col='grey50')
# lines(fit$mu.fv)
#
# fit = gam(y~s(x,bs="ad"),family=nb(link = "log"),data=list(y=y,x=seq(0,1,length.out = n)))
# plot(y,col='grey80')
# lines(exp(b),col='grey50')
# lines(fit$fitted.values)
#
# fit = gam(y~s(x,bs="ad"),family=poisson(link = "log"),data=list(y=y,x=seq(0,1,length.out = n)))
# plot(y,col='grey80')
# lines(exp(b),col='grey50')
# lines(fit$fitted.values)
#
# fit = gam(y~s(x,bs="ad"),family=quasipoisson(link = "log"),data=list(y=y,x=seq(0,1,length.out = n)))
# plot(y,col='grey80')
# lines(exp(b),col='grey50')
# lines(fit$fitted.values)
#
#
# fit = trendfilter(y,family='poisson',k=0)
# plot(y,col='grey80')
# lines(exp(b),col='grey50')
# lines(exp(fit$beta[,15]))
#




