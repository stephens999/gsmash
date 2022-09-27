

library(mr.ash)
poisson_nugget_glm_GMGM = function(X,y,
                                   sa2 = NULL,
                                   beta.init = NULL,
                                   w = NULL,
                                   sigma2 = NULL,
                                   control = list(),
                                   printevery=10){
  n = length(y)
  p = dim(X)[2]
  control <- modifyList(mr_ash_control_default(),control,keep.null = TRUE)

  if(is.null(beta.init)){
    beta = as.vector(double(p))
  }

  if(is.null(sigma2)){
    r = drop(log(y+1) - X %*% beta)
    sigma2 = c(var(r))
  }

  if(is.null(sa2)){
    sa2 = (2^((1:20) / 20) - 1)^2
    #sa2 = sa2 / median(colSums(X^2)) * n
  }
  K = length(sa2)

  if(is.null(w)){
    w = rep(1,K)/K
  }
  Phi = matrix(rep(w,each = p),nrow = p)
  B = matrix(0,nrow=p,ncol=K)
  V = matrix(0.1,nrow=p,ncol=K)
  e = rep(0,n)
  tau2 = rep(0.1,n)

  obj = c(-Inf)

  for(iter in 1:control$max.iter){

    # update b,v
    bv = Update_bv(X,y,B,V,e,tau2,Phi,sa2)
    B = bv$B
    V = bv$V
    # update e, tau2
    etau2 = Update_etau2(X,y,B,V,e,tau2,Phi,sa2)
    e = etau2$e
    tau2 = etau2$tau2

    # update q(z)
    Phi = Update_Phi(X,y,B,V,e,tau2,Phi,w,sigma2,sa2)
    # update sigma2, w
    if(control$update.sigma2){
      sigma2 = mean(e^2+tau2)
    }
    if(control$update.pi){
      w = colMeans(Phi)
      w = pmax(w,1e-15)
    }

    # calc obj
    obj[iter+1] = calc_obj_GMGM(X,y,B,V,e,tau2,Phi,w,sigma2,sa2)
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
  m = drop(rowSums(Phi*B))
  return(list(B=B,V=V,w=w,sigma2=sigma2,Phi=Phi,m=m,obj=obj))

}

Update_Phi  = function(X,y,B,V,e,tau2,Phi,w,sigma2,sa2){
  n = length(y)
  p = dim(X)[2]
  K = length(sa2)
  b = rowSums(B*Phi)
  W = matrix(w,nrow=p,ncol=K,byrow = T)
  Sa2 = matrix(sa2,nrow=p,ncol=K,byrow = T)

  calc_Exb_for_phi = function(b,v,j,X,y,B,V,Phi){
    n = length(y)
    p = dim(X)[2]
    temp = matrix(nrow=n,ncol=p)
    for(i in 1:n){
      temp[i,] = c(rowSums(Phi*exp(X[i,]*B + X[i,]^2*V/2)))
    }
    return(c(apply(temp[,-j],1,prod))*exp(outer(X[,j],b)+outer(X[,j]^2,v)/2))
  }


  temp = matrix(nrow=p,ncol=K)
  for(j in 1:p){
    temp[j,] = c(colSums(y*outer(X[,j],B[j,])-exp(e+tau2/2)*calc_Exb_for_phi(B[j,],V[j,],j,X,y,B,V,Phi)))
  }
  Phi = temp + log(W) - log(Sa2)/2 - (B^2+V)/2/Sa2 + log(V)/2
  Phi = Phi - c(apply(Phi,1,max))
  Phi = exp(Phi)
  Phi = Phi/rowSums(Phi)
  Phi = pmax(Phi,1e-15)
  Phi
}

calc_obj_GMGM = function(X,y,B,V,e,tau2,Phi,w,sigma2,sa2){
  n = length(y)
  K = length(sa2)
  p = dim(X)[2]
  b = rowSums(B*Phi)
  W = matrix(w,nrow=p,ncol=K,byrow = T)
  Sa2 = matrix(sa2,nrow=p,ncol=K,byrow = T)
  obj = sum(y*(X%*%b)+y*e-exp(e+tau2/2)*calc_Exb(X,y,B,V,Phi)) + sum(Phi*(log(W)-log(Sa2)/2-(B^2+V)/2/Sa2)) - n/2*log(sigma2) - sum(e^2+tau2)/2/sigma2 - sum(Phi*log(Phi)) + sum(Phi*log(V))/2 + sum(log(tau2))/2
  return(obj)
}



Update_bv = function(X,y,B,V,e,tau2,Phi,sa2){
  n = length(y)
  p = dim(X)[2]
  K = length(sa2)
  for(j in 1:p){
    for(k in 1:K){
      B[j,k] = opt_b(B[j,k],j,k,X,y,B,V,e,tau2,Phi,sa2)
      V[j,k] = opt_v(V[j,k],j,k,X,y,B,V,e,tau2,Phi,sa2)
    }
  }
  return(list(B=B,V=V))
}

Update_etau2 = function(X,y,B,V,e,tau2,Phi,sa2){
  n = length(y)
  for(i in 1:n){
    e[i] = opt_e(e[i],y[i],X,tau2[i],Phi,B,V,sigma2)
    tau2[i] = opt_tau2(tau2[i],y[i],X,e[i],Phi,B,V,sigma2)
  }
  return(list(e=e,tau2=tau2))
}

opt_tau2 = function(tau2_init,y,X,e,Phi,B,V,sigma2){
  sol = nleqslv::nleqslv(log(tau2_init),opt_tau2_obj,
                         y=y,X=X,
                         e=e,Phi=Phi,B=B,V=V,sigma2=sigma2,method='Newton')$x
  return(exp(sol))
}

opt_tau2_obj = function(tau2_tilde,y,X,e,Phi,B,V,sigma2){
  return(exp(tau2_tilde)*( - exp(e+exp(tau2_tilde)/2)*calc_Exb(X,y,B,V,Phi) - 1/sigma2 + 1/exp(tau2_tilde)))
}


opt_e = function(e_init,y,X,tau2,Phi,B,V,sigma2){
  sol = nleqslv::nleqslv(e_init,opt_e_obj,
                         y=y,X=X,
                         tau2=tau2,Phi=Phi,B=B,V=V,sigma2=sigma2,method='Newton')$x
  return(sol)
}
opt_e_obj = function(e,y,X,tau2,Phi,B,V,sigma2){
  return(y - exp(e+tau2/2)*calc_Exb(X,y,B,V,Phi) - e/sigma2)
}

opt_b = function(b_init,j,k,X,y,B,V,e,tau2,Phi,sa2){
  sol = nleqslv::nleqslv(x=b_init,fn=opt_b_obj,
                         jac=NULL,
                         j=j,k=k,
                         X=X,y=y,
                         B=B,V=V,
                         e=e,tau2=tau2,
                         Phi=Phi,sa2=sa2,method='Newton')$x
  return(sol)
}

opt_v = function(v_init,j,k,X,y,B,V,e,tau2,Phi,sa2){
  sol = nleqslv::nleqslv(x=log(v_init),fn=opt_v_obj,
                         jac=NULL,
                         j=j,k=k,
                         X=X,y=y,
                         B=B,V=V,
                         e=e,tau2=tau2,
                         Phi=Phi,sa2=sa2,method='Newton')$x
  return(exp(sol))
}

calc_Exb = function(X,y,B,V,Phi){
  n = length(y)
  p = dim(X)[2]
  temp = matrix(nrow=n,ncol=p)
  for(i in 1:n){
    temp[i,] = c(rowSums(Phi*exp(X[i,]*B + X[i,]^2*V/2)))
  }
  return(c(apply(temp,1,prod)))

}

calc_Exb_grad = function(b,v,j,X,y,B,V,Phi){
  n = length(y)
  p = dim(X)[2]
  temp = matrix(nrow=n,ncol=p)
  for(i in 1:n){
    temp[i,] = c(rowSums(Phi*exp(X[i,]*B + X[i,]^2*V/2)))
  }
  return(c(apply(temp[,-j],1,prod))*exp(X[,j]*b+X[,j]^2*v/2))
}



opt_b_obj = function(b,j,k,X,y,B,V,e,tau2,Phi,sa2){
  x = X[,j]
  v = V[j,k]
  return(sum(y*x-exp(e+tau2/2)*x*calc_Exb_grad(b,v,j,X,y,B,V,Phi)) - b/sa2[k])
}

opt_v_obj = function(v_tilde,j,k,X,y,B,V,e,tau2,Phi,sa2){
  x = X[,j]
  b = B[j,k]
  return(exp(v_tilde)*(sum(-exp(e+tau2/2)*x^2/2*calc_Exb_grad(b,exp(v_tilde),j,X,y,B,V,Phi)) - 1/2/sa2[k]+1/2/exp(v_tilde)))
}


