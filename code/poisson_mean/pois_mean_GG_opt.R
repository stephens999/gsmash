## previously, I used a for loop to solve for posterior mean, variance for each observation.
## This function is a vectorised version of the for loop + pois_mean_GG1

## This is not faster than the 1-by-1 version. Very likely due to the numerical calc of hessian matrix.
pois_mean_GG_opt = function(x,
                            beta,
                            sigma2,
                            optim_method = 'BFGS',
                            m_init  = NULL,
                            s2_init = NULL){
  n = length(x)
  # init m, s2
  if(is.null(m_init)){
    m = rep(0,n)
  }else{
    m = m_init
  }

  if(is.null(s2_init)){
    s2 = rep(1,n)
  }else{
    s2 = s2_init
  }


  opt = optim(c(m,log(s2)),
              fn = pois_mean_GG_opt_obj,
              gr = pois_mean_GG_opt_obj_gradient,
              x=x,
              beta=beta,
              sigma2=sigma2,
              n=n,
              method = optim_method)

  return(list(m=opt$par[1:n],s2=exp(opt$par[(n+1):(2*n)]),obj=-opt$value))

  # opt = nlm(f_obj_nlm,c(m,log(s2)),x=x,
  #           beta=beta,
  #           sigma2=sigma2,
  #           n=n)
  # return(list(m=opt$estimate[1:n],s2=exp(opt$estimate[(n+1):(2*n)]),obj=-opt$minimum))

}

#'calculate objective function
pois_mean_GG_opt_obj = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+log(exp(v))/2))
}

#'calculate gradient vector
pois_mean_GG_opt_obj_gradient = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}

f_obj_nlm = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  out = -sum(x*m-exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+log(exp(v))/2)
  g1 = -(x-exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  attr(out,'gradient') = c(g1,g2)
  return(out)
}


