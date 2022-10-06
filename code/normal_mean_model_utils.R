library(ashr)

#'log marginal likelihood of normal mean model
#'@param x data vector of length n
#'@param s standard error
#'@param w prior weights
#'@param grid grid of sd in prior
#'@return a vector of length n
l_nm = function(x,s,w,grid){
  return(log(f_nm(x,s,w,grid)))
}

#'@return a n by K matrix of normal density
nm_density = function(x,s,grid){
  K = length(grid)
  n = length(x)
  if(length(s)==1){
    s = rep(s,n)
  }
  sdmat = sqrt(outer(s^2,grid^2,FUN="+"))
  return(stats::dnorm(outer(x,rep(1,K),FUN="*")/sdmat)/sdmat)
}

#'@return a vector of likelihood, of length n
f_nm = function(x,s,w,grid){
  return(c(nm_density(x,s,grid)%*%w))
}

# n = 100000
# x = 1:n
# t1 = Sys.time();a = matrix(x,nrow=n,ncol=200,byrow=FALSE);Sys.time()-t1
# t1 = Sys.time();a = outer(x,rep(1,200));Sys.time()-t1

# l_nm_d1 = function(x,w,grid){
#   s = sqrt(exp(-x))
#   f = sum(w*dnorm(x,0,sqrt(grid^2+s^2)))
#   f_d1 = sum(w*dnorm(x,0,sqrt(grid^2+s^2))*x/(grid^2+s^2))
#   f_d1/f
# }

#'@return a vector of gradient df/dz, length n
f_nm_d1_z = function(x,s,w,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c(-(nm_density(x,s,grid)/vmat*x)%*%w))
}

#'@return a vector of second derivative d^2f/dz^2, length n
f_nm_d2_z = function(x,s,w,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,grid)*(x^2/vmat^2-1/vmat))%*%w))
}

#'@return a vector of third derivative d^3f/dz^3, length n
f_nm_d3_z = function(x,s,w,grid){
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,grid)*(3*x/vmat^2-x^3/vmat^3))%*%w))
}

#'@return a vector of derivative dl_nm/dz, length n
l_nm_d1_z = function(x,s,w,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  return(f_nm_d1_z(x,s,w,grid)/f_nm(x,s,w,grid))
}

#'@return a vector of second derivative d^2l_nm/dz^2, length n
l_nm_d2_z = function(x,s,w,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,grid)
  return(f_nm_d2_z(x,s,w,grid)/temp - (f_nm_d1_z(x,s,w,grid)/temp)^2)
}

#'@return a vector of third derivative d^3l_nm/dz^3, length n
l_nm_d3_z = function(x,s,w,grid){
  if(length(s)==1){
    s = rep(s,length(x))
  }
  temp = f_nm(x,s,w,grid)
  return(f_nm_d3_z(x,s,w,grid)/temp - 3*f_nm_d1_z(x,s,w,grid)*f_nm_d2_z(x,s,w,grid)/temp^2 +2*(f_nm_d1_z(x,s,w,grid)/temp)^3)
}

#'@return a vector of derivative df_nm/ds2, length n
f_nm_d1_s2 = function(x,s,w,grid){
  K = length(grid)
  xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((nm_density(x,s,grid)/vmat^2*(xmat-vmat))%*%w/2))
}

#'@return a vector of second derivative d^2f_nm/dzds2, length n
f_nm_d2_zs2 = function(x,s,w,grid){
  K = length(grid)
  xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  return(c((xmat*(3*vmat-xmat^2)/vmat^3.5*exp(-xmat^2/2/vmat))%*%w/2/sqrt(2*pi)))
}

#'@return a vector of derivative dl_nm/ds2, length n
l_nm_d1_s2 = function(x,s,w,grid){
  return(f_nm_d1_s2(x,s,w,grid)/f_nm(x,s,w,grid))
}

#'@return a vector of second derivative d^2l_nm/dzds2, length n
l_nm_d2_zs2 = function(x,s,w,grid){
  temp = f_nm(x,s,w,grid)
  return(f_nm_d2_zs2(x,s,w,grid)/temp-f_nm_d1_s2(x,s,w,grid)*f_nm_d1_z(x,s,w,grid)/temp^2)
}


#'@return a matrix of gradient, size n* K
f_nm_d1_a = function(x,s,a,grid){
  n = length(x)
  dens_mat = nm_density(x,s,grid)
  return((dens_mat*sum(exp(a)) - c(dens_mat%*%exp(a)))*outer(rep(1,n),exp(a))/sum(exp(a))^2)
}

#'@return a matrix of gradient, size n*K
f_nm_d2_za = function(x,s,a,grid){
  K = length(grid)
  xmat = outer(x,rep(1,K),FUN="*")
  vmat = outer(s^2,grid^2,FUN="+")
  dens_mat = nm_density(x,s,grid)
  lhs = c((dens_mat/vmat)%*%exp(a))
  rhs = (dens_mat/vmat)*sum(exp(a))
  return(outer(x,exp(a))*(lhs-rhs)/sum(exp(a))^2)

}

#'@return a matrix of gradient, size n*K
l_nm_d1_a = function(x,s,a,grid){
  w = softmax(a)
  return(f_nm_d1_a(x,s,a,grid)/f_nm(x,s,w,grid))
}

#'@return a matrix of gradient, size n*K
l_nm_d2_za = function(x,s,a,grid){
  w = softmax(a)
  temp = f_nm(x,s,w,grid)
  return(f_nm_d2_za(x,s,a,grid)/temp - f_nm_d1_a(x,s,a,grid)*f_nm_d1_z(x,s,w,grid)/temp^2)
}




softmax = function(a){
  exp(a-max(a))/sum(exp(a-max(a)))
}

#'posterior mean operator
S = function(x,s,w,grid){
  K = length(w)
  #s = sqrt(exp(-x))
  g = normalmix(pi=w,mean=rep(0,K),sd=grid)
  #return(ashr:::calc_pm())
  fit.ash = ashr::ash(x,s,g=g,fixg=T)
  fit.ash$result$PosteriorMean
}

#'@title inverse operator of S
#'@description  S^{-1}(theta) returns the z such that S(z) = theta
S_inv = function(theta,s,w,grid,z_range=NULL){

  obj = function(z,theta,s,w,grid){
    return(z+s^2*l_nm_d1_z(z,s,w,grid)-theta)
  }
  n = length(theta)
  z_out = double(n)
  for(j in 1:n){
    #print(j)
    if(theta[j]>=0){
      z_out[j] = uniroot(obj,c(theta[j],z_range[2]),
                         theta=theta[j],s=s[j],w=w,grid=grid,extendInt = 'upX')$root
    }else{
      z_out[j] = uniroot(obj,c(z_range[1],theta[j]),
                         theta=theta[j],s=s[j],w=w,grid=grid,extendInt = 'upX')$root
    }

  }
  z_out
}

#'@title compound penalty function of ebnm
nm_penalty_compound = function(theta,s,w,grid){
  return(-l_nm(theta,s,w,grid) - (theta-S(theta,s,w,grid))^2/2/s^2)
}
#'@title gradient of compound penalty function
#'@description -l'(theta) - s^2l'(theta)l''(theta)
nm_penalty_compound_grad = function(theta,s,w,grid){
  return(-l_nm_d1_z(theta,s,w,grid) - s^2*l_nm_d1_z(theta,s,w,grid)*l_nm_d2_z(theta,s,w,grid))
}
#'@title hessian of compound penalty function
#'@description -l''(theta) - s^2(l''(theta)^2+l'''(theta)l'(theta))
nm_penalty_compound_hess = function(theta,s,w,grid){
  return(-l_nm_d2_z(theta,s,w,grid) - s^2*(l_nm_d2_z(theta,s,w,grid)^2+l_nm_d3_z(theta,s,w,grid)*l_nm_d1_z(theta,s,w,grid)))
}

#'@title penalty function of ebnm
nm_penalty = function(theta,s,w,grid,z_range=NULL){
  if(is.null(z_range)){
    z_range = range(theta) + c(-1,1)
  }
  z = S_inv(theta,s,w,grid,z_range)
  original_penalty = nm_penalty_compound(z,s,w,grid)
  return(original_penalty)
}

#'@title gradient of penalty function of ebnm
nm_penalty_grad = function(theta,s,w,grid,z_range=NULL){
  if(is.null(z_range)){
    z_range = range(theta) + c(-1,1)
  }
  z = S_inv(theta,s,w,grid,z_range)
  return((z-theta)/s^2)
}

#'@title hessian of penalty function of ebnm
nm_penalty_hess = function(theta,s,w,grid,z_range=NULL){
  if(is.null(z_range)){
    z_range = range(theta) + c(-1,1)
  }
  z = S_inv(theta,s,w,grid,z_range)
  temp = l_nm_d2_z(z,s,w,grid)
  return(-temp/(1+s^2*temp))
}
