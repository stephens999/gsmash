#'@title Solve Gaussian approximation to Poisson data problem
#'@param x data
#'@param beta prior mean
#'@param sigma2 prior variance
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{m:}{posterior mean}
#'  \item{s2:}{posterior variance}
#'  \item{obj:}{objective function}
#'@examples
#'x = 10
#'beta = 0
#'sigma2 = 100
#'vga_pois(x,beta,sigma2)
#'@details The problem is
#'\deqn{x\sim Poisson(\exp(\mu)),}
#'\deqn{\mu\sim N(\beta,\sigma^2).}
#'@export

vga_pois = function(x,
                    s,
                    beta,
                    sigma2,
                    maxiter=1000,
                    tol=1e-5,
                    m_init  = NULL,
                    s2_init = NULL){
  # init m, s2
  if(is.null(m_init)){
    m_init = ifelse(x==0,0,log(x))
  }
  m = m_init


  if(is.null(s2_init)){
    s2_init = ifelse(x==0,1,1/x)
  }
  s2 = s2_init
  theta = c(m,log(s2))
  obj0 = vga_pois_obj(theta[1],theta[2],x,s,beta,sigma2)
  obj = obj0
  iter = 0
  while(iter<maxiter){
    iter = iter + 1
    #print(theta)
    d = direct_v(theta[1],theta[2],x,s,beta,sigma2)
    #d = d/sqrt(sum(d^2))
    theta = theta - d
    theta[2] = min(theta[2],1)
    obj[iter+1] = vga_pois_obj(theta[1],theta[2],x,s,beta,sigma2)
    if(abs(obj[iter+1]-obj[iter])<tol){
      break
    }
  }
  if(obj[iter+1]<obj0){
    return(list(m=m_init,s2=s2_init,obj=obj0))
  }else{
    return(list(m=theta[1],s2=exp(theta[2]),obj=obj))
  }

}

#'@title calculate Hessian*gradient vector, v = log(s2)
direct_v = function(m,v,x,s,beta,sigma2){
  #browser()
  temp = s*exp(m+exp(v)/2)
  a = -temp-1/sigma2
  b = -temp*exp(v)/2
  c = b
  d = -temp*(exp(2*v)/4+exp(v)/2)-exp(v)/2/sigma2
  g1 = x-temp-m/sigma2+beta/sigma2
  g2 = -temp*exp(v)/2 - exp(v)/2/sigma2 + 1/2
  return(c(d*g1-b*g2,a*g2-c*g1)/(a*d-b*c))
}

#'@title calculate objective function
vga_pois_obj = function(m,v,x,s,beta,sigma2){
  return(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2)
}



#'@title Solve Gaussian approximation to Poisson data problem for posterior mean
#'
vga_pois_mean = function(){

}






