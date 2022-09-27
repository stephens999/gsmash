#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian Mixture prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param beta prior mean
#'@param sigma2 prior variance
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{m:}{posterior mean}
#'  \item{s2:}{posterior variance}
#'  \item{obj:}{objective function values}
#'  \item{beta:}{prior mean}
#'  \item{sigma2:}{prior variance}
#'  @example
#'  n = 10000
#'  mu = rnorm(n)
#'  x = rpois(n,exp(mu))
#'  pois_mean_GMG(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim \sum_k N(\beta,\sigma_k^2).}

pois_mean_GMG = function(x,
                         w = NULL,
                         beta = NULL,
                         sigma2k=NULL,
                         optim_method = 'BFGS',
                         maxiter = 100,
                         tol = 1e-3){
  if(is.null(sigma2k)){
    #sigma2k = c(1e-10,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
    sigma2k = c(1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
  }
  K = length(sigma2k)
  n = length(x)

  if(is.null(beta)){
    beta = log(sum(x)/n)
    est_beta=TRUE
  }else{
    est_beta=FALSE
  }

  m = rep(0,n)
  v = rep(1,n)

  Sigma2k = matrix(sigma2k,nrow=n,ncol=K,byrow=T)
  if(is.null(w)){
    w = rep(1/K,K)
  }

  qz = matrix(0,nrow=n,ncol=K)

  obj = rep(0,maxiter+1)
  obj[1] = -Inf

  for(iter in 1:maxiter){

    qz = matrix(log(w),nrow=n,ncol=K,byrow=T)-log(Sigma2k)/2-matrix(m^2+v-2*m*beta+beta^2,nrow=n,ncol=K,byrow=F)/Sigma2k/2
    qz = qz - apply(qz,1,max)
    qz = exp(qz)
    qz = qz/rowSums(qz)
    qz = pmax(qz,1e-15)

    for(i in 1:n){
      temp = pois_mean_GMG1(x[i],qz[i,],beta,sigma2k,optim_method,m[i],v[i])
      m[i] = temp$m
      v[i] = temp$s2
    }

    if(est_beta){
      beta = sum(qz/Sigma2k*matrix(m,nrow=n,ncol=K,byrow=F)/sum(qz/Sigma2k))
    }

    w = colMeans(qz)
    w = pmax(w, 1e-8)

    #print(w)
    obj[iter+1] = pois_mean_GMG_obj(x,m,v,w,beta,Sigma2k,qz)
    if((obj[iter+1] - obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }


  }

  return(list(m=m,v=v,obj=obj,w=w,qz=qz,beta=beta))

}

pois_mean_GMG_obj = function(x,m,s2,w,beta,Sigma2k,qz){
  n = length(x)
  K = length(w)
  W = matrix(w,nrow=n,ncol=K,byrow=T)
  M = matrix(m,nrow=n,ncol=K,byrow=F)
  S2 = matrix(s2,nrow=n,ncol=K,byrow=F)
  return(sum(x*m-exp(m+s2/2))+sum(qz*(log(W)-log(Sigma2k)/2-(M^2+S2-2*M*beta+beta^2)/2/Sigma2k))-sum(qz*log(qz))+sum(log(s2))/2)
}


pois_mean_GMG1 = function(x,
                          qz,
                          beta,
                          sigma2k,
                          optim_method,
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
              fn = pois_mean_GMG1_obj,
              gr = pois_mean_GMG1_obj_gradient,
              x=x,
              qz=qz,
              beta=beta,
              sigma2k=sigma2k,
              method = optim_method)

  return(list(m=opt$par[1],s2=exp(opt$par[2]),obj=-opt$value))
}


#'calculate objective function
pois_mean_GMG1_obj = function(theta,x,qz,beta,sigma2k){
  return(-(x*theta[1]-exp(theta[1]+exp(theta[2])/2)-(theta[1]^2+exp(theta[2])-2*theta[1]*beta)/2*sum(qz/sigma2k)+log(exp(theta[2]))/2))
}

#'calculate gradient vector
pois_mean_GMG1_obj_gradient = function(theta,x,qz,beta,sigma2k){
  g1 = -(x-exp(theta[1]+exp(theta[2])/2)-theta[1]*sum(qz/sigma2k)+beta*sum(qz/sigma2k))
  g2 = -(-exp(theta[2])/2*exp(theta[1]+exp(theta[2])/2) - exp(theta[2])/2*sum(qz/sigma2k) + 1/2)
  return(c(g1,g2))
}




