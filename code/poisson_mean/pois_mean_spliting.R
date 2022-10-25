
source("code/poisson_mean/pois_mean_GG.R")
pois_mean_split = function(x,s=NULL,sigma2 = NULL,tol=1e-5,maxiter=1e3,ebnm_params=NULL,optim_method ='L-BFGS-B'){
  n = length(x)
  obj = rep(0,maxiter+1)
  obj[1] = -Inf
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }

  b_pm = rep(0,n)
  #b_pv = rep(1/n,n)
  mu_pm = rep(0,n)
  mu_pv = rep(1/n,n)
  if(is.null(sigma2)){
    sigma2 = var(log(x+1))
  }
  for (iter in 1:maxiter) {

    # # VGA
    # for(i in 1:n){
    #   temp = pois_mean_GG1(x[i],s[i],b_pm[i],sigma2,optim_method,mu_pm[i],mu_pv[i])
    #   mu_pm[i] = temp$m
    #   mu_pv[i] = temp$v
    # }

    opt = optim(c(mu_pm,log(mu_pv)),
                fn = pois_mean_GG_opt_obj,
                gr = pois_mean_GG_opt_obj_gradient,
                x=x,
                s=s,
                beta=b_pm,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    mu_pm = opt$par[1:n]
    mu_pv = exp(opt$par[(n+1):(2*n)])

    # EBNM
    res = ebnm(mu_pm,sqrt(sigma2),
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    b_pm = res$posterior$mean
    b_pv = res$posterior$sd^2
    H = res$log_likelihood + n*(log(2*pi*sigma2)/2)+sum((mu_pm^2-2*mu_pm*b_pm+b_pm^2+b_pv)/sigma2/2)

    # Update sigma2
    sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)

    # ELBO
    obj[iter+1] = sum(x*mu_pm-s*exp(mu_pm+mu_pv/2)) - sum(lfactorial(x)) - n/2*log(2*pi*sigma2) - sum(mu_pm^2 + mu_pv + b_pm^2 + b_pv - 2*mu_pm*b_pm)/2/sigma2 + H + sum(log(2*pi*mu_pv))/2 - n/2
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }

  }
  return(list(posteriorMean = mu_pm,posteriorVar = mu_pv,sigma2=sigma2,ebnm_res=res,obj_value=obj))

}


ebnm_params_default = function(){
  return(list(prior_family='point_laplace',
              mode='estimate',
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE,
              output = output_default(),
              optmethod = NULL))
}
