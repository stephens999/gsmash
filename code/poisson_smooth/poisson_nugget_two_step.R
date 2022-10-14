library(funflash)
#'@title Smooth over-dispersed Poisson sequence
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param Eb_init,sigma2_init initial values of smooth mean and nugget effect.
#'@examples
#' set.seed(12345)
#' n=2^9
#' sigma=0.5
#' mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
#' x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
#' fit = vag_pois_smooth(x,maxiter=30)
#' plot(x,col='grey80')
#' lines(exp(fit$Eb))
#' fit$sigma2
#' plot(fit$obj)
#'

vag_pois_smooth = function(x,
                           Eb_init = NULL,
                           sigma2_init = NULL,
                           est_sigma2 = TRUE,
                           maxiter = 100,
                           tol=1e-5,
                           filter.number = 1,
                           family = 'DaubExPhase',
                           printevery = 10){

  n = length(x)
  if(is.null(Eb_init)){
    Eb = log(runmed(x,1 + 2 * min((n-1)%/% 2, ceiling(0.1*n)))+0.01)
  }else{
    Eb = Eb_init
  }
  if(is.null(sigma2_init)){
    sigma2 = var(log(x+0.01)-Eb)
  }else{
    sigma2 = sigma2_init
  }

  m = rep(0,n)
  s2 = rep(1,n)
  obj = -Inf

  for(iter in 1:maxiter){
    # get m, s^2

    # this can be parallel when n is large?
    for(i in 1:n){
      qmu = vag_pois_optim(x[i],Eb[i],sigma2,m_init=m[i],s2_init=s2[i])
      m[i] = qmu$m
      s2[i] = qmu$s2
    }
    # qmu = mclapply(1:n,function(i){
    #   fit = vag_pois(x[i],Eb[i],sigma2,maxiter,tol)
    #   return(c(fit$m,fit$s2))
    # },mc.cores = 4)
    # qmu = do.call(rbind,qmu)
    # m = qmu[,1]
    # s2 = qmu[,2]
    # get Eb, Eb2
    qb = smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family)
    Eb = qb$mu.est
    Eb2 = qb$mu.est.var + Eb^2
    # get sigma2
    if(est_sigma2){
      sigma2 = mean(m^2+s2+Eb2-2*m*Eb)
    }


    # calc obj
    obj[iter+1] = vag_pois_smooth_obj(x,m,s2,Eb,Eb2,sigma2,qb$dKL)

    if(iter%%printevery==0){
      print(paste("Done iter",iter,"obj =",obj[iter+1]))
    }

    if((obj[iter+1]-obj[iter])<tol){
      break
    }

  }
  return(list(m=m,s2=s2,Eb=Eb,Eb2=Eb2,sigma2=sigma2,obj=obj))
}

vag_pois_smooth_obj = function(x,m,s2,Eb,Eb2,sigma2,KLb){
  return(sum(x*m-exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb)
}



#'@title Solve Gaussian approximation to Poisson data problem, using optim function for optimization
#'@param x data
#'@param beta prior mean
#'@param sigma2 prior variance
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{m:}{posterior mean}
#'  \item{s2:}{posterior variance}
#'  \item{obj:}{objective function}
#'  @example
#'  x = 10
#'  beta = 0
#'  sigma2 = 100
#'  vag_pois_optim(x,beta,sigma2)
#'@details The problem is
#'\deqn{x\sim Poisson(\exp(\mu)),}
#'\deqn{\mu\sim N(\beta,\sigma^2).}


vag_pois_optim = function(x,
                          beta,
                          sigma2,
                          optim_method = 'BFGS',
                          m_init  = NULL,
                          s2_init = NULL){
  # init m, s2
  if(is.null(m_init)){
    m = 0
  }else{
    m = m_init
  }

  if(is.null(s2_init)){
    s2 = 1
  }else{
    s2 = s2_init
  }

  opt = optim(c(m,log(s2)),
              fn = vag_pois_obj_optim,
              gr = gradient_optim,
              x=x,
              beta=beta,
              sigma2=sigma2,
              method = optim_method)

  return(list(m=opt$par[1],s2=exp(opt$par[2]),obj=-opt$value))
}

#'calculate objective function
vag_pois_obj_optim = function(theta,x,beta,sigma2){
  return(-(x*theta[1]-exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2/sigma2+log(exp(theta[2]))/2))
}

#'calculate gradient vector
gradient_optim = function(theta,x,beta,sigma2){
  g1 = -(x-exp(theta[1]+exp(theta[2])/2)-theta[1]/sigma2+beta/sigma2)
  g2 = -(-exp(theta[2])/2*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2/sigma2 + 1/2)
  return(c(g1,g2))
}

