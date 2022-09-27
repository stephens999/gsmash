
library(mr.ash)
poisson_nugget_glm_GMG = function(X,y,
                                   sa2 = NULL,
                                   beta.init = NULL,
                                   w = NULL,
                                   sigma2 = NULL,
                                   control = list(),
                                  printevery = 10){

  # init things
  # to be done


  n = length(y)
  p = dim(X)[2]
  if(is.null(sa2)){
    sa2 = (2^((1:20) / 20) - 1)^2
    #sa2 = sa2 / median(colSums(X^2)) * n
  }
  K = length(sa2)

  control = modifyList(mr_ash_control_default(),control,keep.null = TRUE)

  if(is.null(beta.init)){
    m = as.vector(double(p))
  }else{
    m = beta.init
  }

  if(is.null(sigma2)){
    r = drop(log(y+1) - X %*% m)
    sigma2 = c(var(r))
  }

  if(is.null(w)){
    w = rep(1,K)/K
  }

  Phi = matrix(rep(w,each = p),nrow = p)
  Sa2 = matrix(sa2,nrow=p,ncol=K,byrow=TRUE)
  X2 = X^2
  obj = c(-Inf)

  v = rep(0.1,p)
  e = rep(0,n)
  d = rep(1,n)

  for(iter in 1:control$max.iter){
    # update m, v
    m = update_m(X,y,v,e,d,Phi,Sa2,X2,m)
    v = update_v(X,y,m,e,d,Phi,Sa2,X2,v)
    # update e,d
    e = update_e(X,y,m,v,d,sigma2,X2,e)
    d = update_d(X,y,m,v,e,sigma2,X2,d)
    # update gamma
    W = matrix(w,nrow=p,ncol=K,byrow=TRUE)
    Phi = log(W) - matrix(m^2+v,nrow=p,ncol=K,byrow=FALSE)/2/Sa2
    Phi = Phi - c(apply(Phi,1,max))
    Phi = exp(Phi)
    Phi = Phi/c(rowSums(Phi))
    Phi = pmax(Phi,1e-15)
    # update w,sigma2
    if(control$update.sigma2){
      sigma2 = mean(e^2+d)
    }
    if(control$update.pi){
      w = colMeans(Phi)
      w = pmax(w,1e-15)
    }


    # calc obj
    obj[iter+1] = calc_obj_GMG(X,y,m,v,e,d,Phi,w,sigma2,Sa2,X2)
    if(iter%%printevery==0){
      print(sprintf("At iter %1.0f, ELBO=%.3f",iter,obj[iter+1]))
    }
    if((obj[iter+1]-obj[iter])<0){
      message('An iteration decreases ELBO')
    }
    if(abs(obj[iter+1]-obj[iter])<control$convtol){
      break
    }
  }
  return(list(m=m,v=v,e=e,d=d,w=w,sigma2=sigma2,obj=obj,Phi=Phi))

}

calc_obj_GMG = function(X,y,m,v,e,d,Phi,w,sigma2,Sa2,X2){
  n = length(y)
  p = dim(X)[2]
  K = length(w)
  W = matrix(w,nrow=p,ncol=K,byrow=TRUE)
  a = get_a(Phi,Sa2)
  val = t(y)%*%X%*%m + t(y)%*%e - sum(exp(X%*%m+X2%*%v/2+e+d/2)) + sum(Phi*(log(W)-log(Sa2)/2)) - sum(a*m^2)/2-sum(a*v)/2 - n*log(sigma2)/2 - (sum(e^2)+sum(d))/2/sigma2 + sum(log(v))/2 + sum(log(d))/2 - sum(Phi*log(Phi))
  return(drop(val))
}

update_d = function(X,y,m,v,e,sigma2,X2,d_init,fixed_point_iter = 1000,tol=1e-12){
  d_old = d_init
  for(fiter in 1:fixed_point_iter){
    d_new = 1/(1/sigma2+exp(X%*%m+X2%*%v/2+e+d_old/2))
    if(sqrt(sum((d_new-d_old)^2))<tol){
      break
    }else{
      d_old = d_new
    }
  }
  return(drop(d_new))
}

G_d = function(d,X,X2,m,v,e,sigma2){
  return(c(-1/2*exp(X%*%m + X2%*%v/2+e+d/2) - 1/2/sigma2+1/2/d))
}

G_v = function(v,X,X2,m,e,d,Phi,Sa2){
  a = get_a(Phi,Sa2)
  return(c(-1/2*t(X2)%*%exp(X%*%m+X2%*%v/2+e+d/2)-a/2+1/2/v))

}

update_e = function(X,y,m,v,d,sigma2,X2,e_init){
  sol = nleqslv::nleqslv(e_init,G_e,H_e,
                         X=X,y=y,m=m,v=v,d=d,sigma2=sigma2,X2=X2,method='Newton')$x
  return(sol)
}

G_e = function(e,X,y,m,v,d,sigma2,X2){
  grad = y - drop(exp(X%*%m + X2%*%v/2+e+d/2)) - e/sigma2
  return(grad)
}

H_e = function(e,X,y,m,v,d,sigma2,X2){
  hess = -diag(drop(exp(X%*%m+X2%*%v/2+e+d/2)) + 1/sigma2)
  return(hess)
}

update_v = function(X,y,m,e,d,Phi,Sa2,X2,v_init,fixed_point_iter = 1000,tol=1e-12){
  v_old = v_init
  a = get_a(Phi,Sa2)
  for(fiter in 1:fixed_point_iter){
    v_new = 1/(a+t(X2)%*%exp(X%*%m+X2%*%v_old/2+e+d/2))
    if(sqrt(sum((v_new-v_old)^2))<tol){
      break
    }else{
      v_old = v_new
    }
  }
  return(drop(v_new))
}

update_m = function(X,y,v,e,d,Phi,Sa2,X2,m_init){
  sol = nleqslv::nleqslv(m_init,G_m,H_m,
                         X=X,y=y,v=v,e=e,d=d,Phi=Phi,Sa2=Sa2,X2=X2,method='Newton')$x
  return(sol)
}

get_a = function(Phi,Sa2){
  return(c(rowSums(Phi/Sa2)))
}

G_m = function(m,X,y,v,e,d,Phi,Sa2,X2){
  a = get_a(Phi,Sa2)
  grad = t(X)%*%y - t(X)%*%exp(X%*%m + X2%*%v/2 + e + d/2) - a*m
  return(grad)
}

H_m = function(m,X,y,v,e,d,Phi,Sa2,X2){
  a = get_a(Phi,Sa2)
  #hess = -t(diag(c(exp(X%*%m + X2%*%v/2+e+d/2)))%*%X)%*%X - diag(a)
  hess = -t(c(exp(X%*%m + X2%*%v/2+e+d/2))*X)%*%X - diag(a)
  return(hess)
}
