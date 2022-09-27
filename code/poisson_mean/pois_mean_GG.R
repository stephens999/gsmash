#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param beta prior mean
#'@param sigma2 prior variance
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{m:}{posterior mean}
#'  \item{v:}{posterior variance}
#'  \item{obj:}{objective function values}
#'  \item{beta:}{prior mean}
#'  \item{sigma2:}{prior variance}
#'  @example
#'  n = 10000
#'  mu = rnorm(n)
#'  x = rpois(n,exp(mu))
#'  pois_mean_GG(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}

pois_mean_GG = function(x,
                        beta = NULL,
                        sigma2=NULL,
                        optim_method = 'BFGS',
                        maxiter = 100,
                        tol = 1e-3){

  # init the posterior mean and variance?
  n = length(x)
  m = log(x+0.1)
  v = rep(1,n)

  #
  if(is.null(beta) | is.null(sigma2)){

    if(is.null(beta)){
      est_beta = TRUE
    }else{
      est_beta = FALSE
    }
    if(is.null(sigma2)){
      est_sigma2=TRUE
    }else{
      est_sigma2 = FALSE
    }

    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    for(iter in 1:maxiter){
      if(est_beta){
        beta = mean(m)
      }
      if(est_sigma2){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
      }
      for(i in 1:n){
        temp = pois_mean_GG1(x[i],beta,sigma2,optim_method,m[i],v[i])
        m[i] = temp$m
        v[i] = temp$v
      }
      obj[iter+1] = pois_mean_GG_obj(x,beta,sigma2,m,v)
      if((obj[iter+1] - obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }

  }else{

    for(i in 1:n){
      temp = pois_mean_GG1(x[i],beta,sigma2,optim_method,m[i],v[i])
      m[i] = temp$m
      v[i] = temp$v
    }

    obj = pois_mean_GG_obj(x,beta,sigma2,m,v)

  }

  return(list(m=m,beta=beta,sigma2=sigma2,v=v,obj=obj))

}


pois_mean_GG_obj = function(x,beta,sigma2,m,v){
  return(sum(x*m-exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2))
}

#'@param x a data point
#'@param beta prior mean
#'@param sigma2 prior variance
#'@param optim_method optimization method in `optim` function
pois_mean_GG1 = function(x,
                            beta,
                            sigma2,
                            optim_method = 'BFGS',
                            m_init  = NULL,
                            v_init = NULL){
  # init m, v
  if(is.null(m_init)){
    m = 0
  }else{
    m = m_init
  }

  if(is.null(v_init)){
    v = 1
  }else{
    v = v_init
  }

  opt = optim(c(m,log(v)),
              fn = pois_mean_GG1_obj,
              gr = pois_mean_GG1_obj_gradient,
              x=x,
              beta=beta,
              sigma2=sigma2,
              method = optim_method)

  return(list(m=opt$par[1],v=exp(opt$par[2]),obj=-opt$value))
}

#'calculate objective function
pois_mean_GG1_obj = function(theta,x,beta,sigma2){
  return(-(x*theta[1]-exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2/sigma2+log(exp(theta[2]))/2))
}

#'calculate gradient vector
pois_mean_GG1_obj_gradient = function(theta,x,beta,sigma2){
  g1 = -(x-exp(theta[1]+exp(theta[2])/2)-theta[1]/sigma2+beta/sigma2)
  g2 = -(-exp(theta[2])/2*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2/sigma2 + 1/2)
  return(c(g1,g2))
}



