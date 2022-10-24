#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Poisson mean problem, log(1+exp(x)) link function.
#'@param x data vector
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{m:}{posterior mean}
#'  \item{v:}{posterior variance}
#'  \item{obj:}{objective function values}
#'  \item{ebnm_res:}{fitted object from `ebnm`}
#'  @example
#'  n = 10000
#'  mu = rnorm(n)
#'  x = rpois(n,exp(mu))
#'  pois_mean_log1exp(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\log(1+\exp(\mu_i))),}
#'\deqn{\mu_i\sim g(\cdot).}
library(ebnm)
pois_mean_log1exp = function(x,ebnm_params = NULL,tol=1e-5,maxiter=1e3){
  n = length(x)
  kapa = 0.25+0.17*max(x)
  pseudo_s = sqrt(1/kapa)

  mu_tilde = 0
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
  for(iter in 1:maxiter){
    # update g,q by performing ebnm
    pseudo_x = mu_tilde-nll_d1(x,mu_tilde)/kapa
    res = ebnm(pseudo_x,pseudo_s,
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    m = res$posterior$mean
    v = res$posterior$sd^2
    # get Elog(g/q)
    H = res$log_likelihood + sum(log(2*pi*pseudo_s^2)/2)+sum((pseudo_x^2-2*pseudo_x*m+m^2+v)/pseudo_s^2/2)

    # calc objective function
    obj[iter+1] = -sum(nll(x,mu_tilde)+nll_d1(x,mu_tilde)*(m-mu_tilde)+kapa/2*(m^2+v+mu_tilde^2-2*m*mu_tilde))+H
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }
    # update mu_tilde
    mu_tilde = m
  }
  return(list(posteriorMean=m,posteriorVar=v,obj_value=obj,ebnm_res=res,kappa=kapa))
}

nll = function(x,mu){
  return(log(1+exp(mu))-x*log(log(1+exp(mu))))
}

nll_d1 = function(x,mu){
  n = exp(mu)*(log(1+exp(mu))-x)
  d = (1+exp(mu))*log((1+exp(mu)))
  return(n/d)
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



