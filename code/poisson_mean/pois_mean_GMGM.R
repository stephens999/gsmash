#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian Mixture prior, Gaussian Mixture posterior in Poisson mean problem.
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
#'\deqn{\mu_i\sim \sum_k N(\beta,\sigma_k^2).}

source("code/poisson_mean/pois_mean_GG.R")

pois_mean_GMGM = function(x,
                          w = NULL,
                          beta = NULL,
                          sigma2k=NULL,
                          optim_method = 'BFGS',
                          maxiter = 100,
                          tol = 1e-3){
  if(is.null(sigma2k)){

    ## how to choose grid in this case?

    #sigma2k = c(1e-10,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
    sigma2k = c(1e-3, 1e-2, 1e-1, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
  }
  K = length(sigma2k)
  n = length(x)

  # init prior mean,
  if(is.null(beta)){
    beta = log(sum(x)/n)
    est_beta = TRUE
  }else{
    est_beta = FALSE
  }


  M = matrix(0,nrow=n,ncol=K)
  V = matrix(1,nrow=n,ncol=K)
  Sigma2k = matrix(sigma2k,nrow=n,ncol=K,byrow=T)
  X = matrix(x,nrow=n,ncol=K,byrow=F)

  # init prior weights
  if(is.null(w)){
    w = rep(1/K,K)
  }

  qz = matrix(0,nrow=n,ncol=K)

  obj = rep(0,maxiter+1)
  obj[1] = -Inf

  for(iter in 1:maxiter){

    # update posterior mean, variances

    ## this can be paralleled?
    for(i in 1:n){
      for (k in 1:K) {
        temp = pois_mean_GG1(x[i],beta,sigma2k[k],optim_method,M[i,k],V[i,k])
        M[i,k] = temp$m
        V[i,k] = temp$v
      }
    }

    # update posterior weights

    qz = X*M-exp(M+V/2)+matrix(log(w),nrow=n,ncol=K,byrow=T)-log(Sigma2k)/2-(M^2+V-2*M*beta+beta^2)/Sigma2k/2 + log(V)/2
    qz = qz - apply(qz,1,max)
    qz = exp(qz)
    qz = qz/rowSums(qz)
    qz = pmax(qz,1e-15)

    # update prior

    if(est_beta){
      beta = sum(M*qz/Sigma2k)/sum(qz/Sigma2k)
    }

    w = colMeans(qz)
    w = pmax(w, 1e-8)


    obj[iter+1] = pois_mean_GMGM_obj(X,M,V,w,beta,Sigma2k,qz)
    if((obj[iter+1] - obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }

  }

  return(list(m=rowSums(qz*M),M=M,V=V,obj=obj,w=w,qz=qz,beta=beta))

}

pois_mean_GMGM_obj = function(X,M,V,w,beta,Sigma2k,qz){
  n = dim(X)[1]
  K = length(w)
  lW = matrix(log(w),nrow=n,ncol=K,byrow=T)
  return(sum(qz*(X*M-exp(M+V/2)+lW-log(Sigma2k)/2-(M^2+V-2*M*beta+beta^2)/2/Sigma2k-log(qz)+log(V)/2)))
}






