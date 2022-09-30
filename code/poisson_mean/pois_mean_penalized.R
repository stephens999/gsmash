source("code/normal_mean_model_utils.R")
#'@title Poisson mean problem via penalized form

pois_mean_penalized_optim = function(x,
                               w=NULL,
                               grid,
                               est_w = TRUE,
                               z_init = NULL,
                               lambda_init = NULL,
                               opt_method='BFGS'){
  n = length(x)
  K = length(grid)
  if(is.null(w)){
    w = rep(1/K,K)
    est_w = TRUE
  }
  if(is.null(z_init)){
    z_init = log(1+x)
  }
  if(is.null(lambda_init)){
    lambda_init = rep(1,n)
  }
  if(est_w){
    fit = optim(c(z_init,-z_init,lambda_init,w),
                L,
                L_grad,
                method=opt_method,
                y=x,
                grid=grid)
    z=fit$par[1:n]
    s2=exp(fit$par[(n+1):(2*n)])
    w = softmax(fit$par[-(1:(3*n))])
    return(list(z=z,
                m=S(z,sqrt(s2),w,grid),
                s2=s2,
                w=w,
                lambda = fit$par[(2*n+1):(3*n)],
                fit=fit))
  }else{
    fit = optim(c(z_init,-z_init,lambda_init),
                L_known_g,
                L_grad_known_g,
                method=opt_method,
                y=x,
                w=w,
                grid=grid)
    z=fit$par[1:n]
    s2=exp(fit$par[(n+1):(2*n)])
    return(list(z=z,
                m=S(z,sqrt(s2),w,grid),
                s2=s2,
                lambda = fit$par[(2*n+1):(3*n)],
                w = w,
                fit=fit))
  }

}

pois_mean_penalized_nleqslv = function(x,
                                     w=NULL,
                                     grid,
                                     est_w = TRUE,
                                     z_init = NULL,
                                     lambda_init = NULL,
                                     opt_method='BFGS'){
  n = length(x)
  K = length(grid)
  if(is.null(w)){
    w = rep(1/K,K)
    est_w = TRUE
  }
  if(is.null(z_init)){
    z_init = log(1+x)
  }
  if(is.null(lambda_init)){
    lambda_init = rep(1,n)
  }
  if(est_w){
    theta0 = c(z_init,-z_init,lambda_init,w)
    fit  = nleqslv(theta0,
                    L_grad,
                    y=x,
                    grid=grid,
                   method = opt_method)
    z=fit$x[1:n]
    s2=exp(fit$x[(n+1):(2*n)])
    w = softmax(fit$x[-(1:(3*n))])
    return(list(z=z,
                m=S(z,sqrt(s2),w,grid),
                s2=s2,
                w=w,
                lambda = fit$x[(2*n+1):(3*n)],
                fit=fit))
  }else{
    theta0 = c(z_init,-z_init,lambda_init)
    fit  = nleqslv(theta0,
                   L_grad_known_g,
                   y=x,
                   w=w,
                   grid=grid,
                   method = opt_method)
    z=fit$x[1:n]
    s2=exp(fit$x[(n+1):(2*n)])
    return(list(z=z,
                m=S(z,sqrt(s2),w,grid),
                s2=s2,
                w=w,
                lambda = fit$x[(2*n+1):(3*n)],
                fit=fit))
  }

}

#'@title Lagrangian of possion mean penalized problem
#'@param theta = (z,v=log(s2),lambda,a)
#'@param y data vector
#'@param grid prior sd
L = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = exp(theta[(n+1):(2*n)])
  lambda = theta[(2*n+1):(3*n)]
  a = theta[-(1:(3*n))]
  return(h_obj_calc(z,s2,a,y,grid) - sum(lambda*h_cstr_calc(z,s2,a,grid)))
}

L_grad = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  v = theta[(n+1):(2*n)]
  s2 = exp(v)
  lambda = theta[(2*n+1):(3*n)]
  a = theta[-(1:(3*n))]
  L_dz = h_obj_d1_z_calc(z,s2,a,y,grid) - lambda*h_cstr_d1_z_calc(z,s2,a,grid)
  L_dv = (h_obj_d1_s2_calc(z,s2,a,y,grid) - lambda*h_cstr_d1_s2_calc(z,s2,a,grid))*exp(v)
  L_da = h_obj_d1_a_calc(z,s2,a,y,grid) - colSums(lambda*h_cstr_d1_a_calc(z,s2,a,grid))
  L_dlambda = h_cstr_calc(z,s2,a,grid)
  return(c(L_dz,L_dv,L_dlambda,L_da))
}

#'@title Lagrangian of possion mean penalized problem, when g is known
# theta = (z,v=log(s2),lambda)
L_known_g = function(theta,y,w,grid){
  n = length(y)
  z = theta[1:n]
  s2 = exp(theta[(n+1):(2*n)])
  lambda = theta[(2*n+1):(3*n)]
  return(h_obj_calc(z,s2,w,y,grid) - sum(lambda*h_cstr_calc(z,s2,w,grid)))
}

L_grad_known_g = function(theta,y,w,grid){
  n = length(y)
  z = theta[1:n]
  v = theta[(n+1):(2*n)]
  s2 = exp(v)
  lambda = theta[(2*n+1):(3*n)]
  a = w
  #browser()
  L_dz = h_obj_d1_z_calc(z,s2,a,y,grid) - lambda*h_cstr_d1_z_calc(z,s2,a,grid)
  L_dv = (h_obj_d1_s2_calc(z,s2,a,y,grid) - lambda*h_cstr_d1_s2_calc(z,s2,a,grid))*exp(v)
  L_dlambda = h_cstr_calc(z,s2,a,grid)
  return(c(L_dz,L_dv,L_dlambda))
}

#'objective function, sum over h_i()
#'@param theta (z,s2,a)
#'@param y data vector
#'@param grid prior sds
#'@return objevtive value, a scalar
h_obj = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_obj_calc(z,s2,a,y,grid))
}

h_obj_calc = function(z,s2,a,y,grid){
  w = softmax(a)
  l_dz = l_nm_d1_z(z,sqrt(s2),w,grid)
  return(sum(exp(z+s2*l_dz)-(y-0.5)*(z+s2*l_dz)-l_nm(z,sqrt(s2),w,grid)-s2*l_dz^2/2))
}

#'objective function derivative wrt z
#'@param theta (z,s2,a)
#'@param y data vector
#'@param grid prior sds
#'@return a length n vector of gradients
h_obj_d1_z = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_obj_d1_z_calc(z,s2,a,y,grid))
}

h_obj_d1_z_calc = function(z,s2,a,y,grid){
  w = softmax(a)
  l_dz = l_nm_d1_z(z,sqrt(s2),w,grid)
  l_dz2 = l_nm_d2_z(z,sqrt(s2),w,grid)
  return(exp(z+s2*l_dz)*(1+s2*l_dz2)-(y-0.5)*(1+s2*l_dz2)-l_dz-s2*l_dz*l_dz2)
}

#'objective function derivative wrt s2
#'@param theta (z,s2,a)
#'@param y data vector
#'@param grid prior sds
#'@return a length n vector of gradients
h_obj_d1_s2 = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_obj_d1_s2_calc(z,s2,a,y,grid))

}

h_obj_d1_s2_calc = function(z,s2,a,y,grid){
  w = softmax(a)
  l_dz = l_nm_d1_z(z,sqrt(s2),w,grid)
  l_dzds2 = l_nm_d2_zs2(z,sqrt(s2),w,grid)
  l_ds2 = l_nm_d1_s2(z,sqrt(s2),w,grid)
  return((exp(z+s2*l_dz)-y+0.5)*(l_dz+s2*l_dzds2)-l_ds2-l_dz^2/2-s2*l_dz*l_dzds2)
}

#'objective function derivative wrt a
#'@param theta (z,s2,a)
#'@param y data vector
#'@param grid prior sds
#'@return a length K vector of gradients
h_obj_d1_a = function(theta,y,grid){
  n = length(y)
  #K = length(grid)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_obj_d1_a_calc(z,s2,a,y,grid))
}

h_obj_d1_a_calc = function(z,s2,a,y,grid){
  w = softmax(a)
  K = length(grid)
  l_dz = l_nm_d1_z(z,sqrt(s2),w,grid)
  l_dadz = l_nm_d2_za(z,sqrt(s2),a,grid)
  l_da = l_nm_d1_a(z,sqrt(s2),a,grid)
  return(colSums((exp(z+s2*l_dz)-y+0.5-l_dz)*s2*l_dadz-l_da))
}

#'@title calculate constraint function
#'@param theta (z,s2,a)
#'@return a vector of length n, constraint function values
h_cstr = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_cstr_calc(z,s2,a,grid))
}

h_cstr_calc = function(z,s2,a,grid){
  w = softmax(a)
  l_dz = l_nm_d1_z(z,sqrt(s2),w,grid)
  return(log(s2)+z+s2*l_dz)
}

#'@title calculate constraint function derivative wrt z
#'@param theta (z,s2,a)
#'@return a vector of length n gradient
h_cstr_d1_z = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_cstr_d1_z_calc(z,s2,a,grid))
}

h_cstr_d1_z_calc = function(z,s2,a,grid){
  w = softmax(a)
  return(1+s2*l_nm_d2_z(z,sqrt(s2),w,grid))
}

#'@title calculate constraint function derivative wrt s2
#'@param theta (z,s2,a)
#'@return a vector of length n gradient
h_cstr_d1_s2 = function(theta,y,grid){
  n = length(y)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_cstr_d1_s2_calc(z,s2,a,grid))
}

h_cstr_d1_s2_calc = function(z,s2,a,grid){
  w = softmax(a)
  return(1/s2+l_nm_d1_z(z,sqrt(s2),w,grid)+s2*l_nm_d2_zs2(z,sqrt(s2),w,grid))
}


#'@title The gradient of constraint function w.r.t a
#'@return a n*K matrix
h_cstr_d1_a = function(theta,y,grid){
  n = length(y)
  K = length(grid)
  z = theta[1:n]
  s2 = theta[(n+1):(2*n)]
  a = theta[-(1:(2*n))]
  return(h_cstr_d1_a_calc(z,s2,a,grid))
}

h_cstr_d1_a_calc = function(z,s2,a,grid){
  w = softmax(a)
  l_dzda = l_nm_d2_za(z,sqrt(s2),a,grid)
  return(s2*l_dzda)
}





#' #'objective function
#' #'@param theta (z,s2,w)
#' #'@param x data vector
#' #'@param grid prior sds
#' f_obj_cstrOpt = function(theta,x,grid){
#'   n = length(x)
#'   z = theta[1:n]
#'   s2 = theta[(n+1):(2*n)]
#'   w = softmax(theta[-(1:(2*n))])
#'   l_nm_d1_z = l_nm_d1(z,sqrt(s2),w,grid)
#'   return(sum(exp(z+s2*l_nm_d1_z)-(x-0.5)*(z+s2*l_nm_d1_z)-l_nm(z,sqrt(s2),w,grid)-s2*l_nm_d1_z^2/2))
#' }
#'
#'
#' f_cstr_cstrOpt = function(theta,x,grid){
#'   n = length(x)
#'   z = theta[1:n]
#'   s2 = theta[(n+1):(2*n)]
#'   w = softmax(theta[-(1:(2*n))])
#'   l_nm_d1_z = l_nm_d1(z,sqrt(s2),w,grid)
#'   return(s2-exp(-z-s2*l_nm_d1_z))
#' }


#' #'objective function
#' #'@param theta (z,w)
#' #'@param x data vector
#' #'@param grid prior sds
#' f_obj_uncstrOpt = function(theta,x,grid){
#'   n = length(x)
#'   z = theta[1:n]
#'   w = softmax(theta[-(1:n)])
#'   s2 = solve_for_s2_given_zg(z,w,grid)
#'   l_nm_d1_z = l_nm_d1(z,sqrt(s2),w,grid)
#'   return(sum(exp(z+s2*l_nm_d1_z)-(x-0.5)*(z+s2*l_nm_d1_z)-l_nm(z,sqrt(s2),w,grid)-s2*l_nm_d1_z^2/2))
#'
#' }

solve_for_s2_given_zg = function(z,w,grid){
  f = function(theta,z,w,grid){
    return(theta+z+exp(theta)*l_nm_d1(z,sqrt(exp(theta)),w,grid))
  }
  n = length(z)
  s2 = rep(0,n)
  for(i in 1:n){
    s2[i] = uniroot(f,c(1e-12,exp(-z[i])),z=z[i],w=w,grid=grid)$root
  }
  return(exp(s2))
}

solve_for_z_given_s2g = function(s2,w,grid){
  f =  function(z,s2,w,grid){
    return(z+log(s2)+s2*l_nm_d1_z(z,sqrt(s2),w,grid))
  }
  n = length(s2)
  z = rep(0,n)
  for(i in 1:n){
    z[i] = uniroot(f,c(-5,5),s2=s2[i],w=w,grid=grid)$root
  }
  return(z)
}



# f_obj_mu = function(mu,y,w,grid){
#   s = sqrt(exp(-mu))
#   l_nm_d1_mu = l_nm_d1(mu,s,w,grid)
#   return(sum(exp(mu+s^2*l_nm_d1_mu)-(y-0.5)*(mu+s^2*l_nm_d1_mu)-l_nm(mu,s,w,grid)-s^2*l_nm_d1_mu^2/2))
#
# }





