
## examples
library(vebpm)
library(fastGHQuad)
library(ebnm)
library(Rfast)
binomial_mean_splitting = function(x,nb=1,
                                   tol=1e-5,
                                   maxiter=1e3,
                                   ebnm_params=NULL,
                                   optim_method ='L-BFGS-B',
                                   n_gh = 20,
                                   printevery = 100,
                                   b_pm_init = NULL,
                                   sigma2_init = NULL,
                                   convergence_criteria = 'obj'){
  n = length(x)
  if(length(nb)==1){
    nb = rep(nb,n)
  }
  obj = rep(0,maxiter+1)
  obj[1] = -Inf
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }

  # const in objective function
  const = sum(lfactorial(nb)-lfactorial(x)-lfactorial(nb-x))

  gh_points = gaussHermiteData(n_gh)

  if(is.null(b_pm_init)){
    b_pm = rep(0,n)
  }else{
    b_pm = b_pm_init
  }

  #b_pv = rep(1/n,n)
  mu_pm = logit(x/nb)
  mu_pm[mu_pm==-Inf] = logit(0.1)
  mu_pm[mu_pm==Inf] = logit(0.9)
  mu_pv = rep(1/n,n)
  if(is.null(sigma2_init)){
    sigma2 = var(mu_pm-b_pm)
  }else{
    sigma2 = sigma2_init
  }
  t_start = Sys.time()
  for (iter in 1:maxiter) {

    opt = vga_binomial(c(mu_pm,log(mu_pv)),x,nb,b_pm,sigma2,gh_points=gh_points)
    mu_pm = opt$m
    mu_pv = opt$v
    # EBNM
    res = ebnm(mu_pm,sqrt(sigma2),
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    b_pm = res$posterior$mean
    b_pv = res$posterior$sd^2
    H = res$log_likelihood + n*(log(2*pi*sigma2)/2)+sum((mu_pm^2-2*mu_pm*b_pm+b_pm^2+b_pv)/sigma2/2)

    # Update sigma2
    sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)

    # ELBO
    obj[iter+1] = binom_mean_split_obj(x,nb,mu_pm,mu_pv,gh_points,const,sigma2,b_pm,b_pv,H)

    if(iter%%printevery==0){
      print(paste('At iter', iter, 'elbo=',round(obj[iter+1],3),'sigma2=',round(sigma2,3)))
    }
    if(convergence_criteria=='obj'){
      if((obj[iter+1]-obj[iter])<tol){
        obj = obj[1:(iter+1)]
        if((obj[iter+1]-obj[iter])<0){
          warning('An iteration decreases ELBO. This is likely due to numerical issues.')
        }
        break
      }
    }
    if(convergence_criteria=='objabs'){
      if(abs(obj[iter+1]-obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }


  }
  t_end = Sys.time()
  return(list(posterior = list(mean_logit = mu_pm,
                               mean_b = b_pm,
                               mean = Esigmoid(mu_pm,mu_pv,gh_points)),
              fitted_g = list(sigma2=sigma2,g_b = res$fitted_g),
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = res,
              run_time = difftime(t_end,t_start,units='secs')))

}


binom_mean_split_obj = function(x,nb,mu_pm,mu_pv,gh_points,const,sigma2,b_pm,b_pv,H){
  n = length(x)
  term1 =sum(x*mu_pm-nb*Elog1pexp(mu_pm,mu_pv,gh_points)- (mu_pm^2 + mu_pv + b_pm^2 + b_pv - 2*mu_pm*b_pm)/2/sigma2+ (log(2*pi*mu_pv))/2)
  term2 =const - n/2*log(2*pi*sigma2) + H  - n/2
  return(term1 + term2)
}

vga_binomial = function(init_val,x,nb,beta,sigma2,method='lbfgs',gh_points){
  if(method!="lbfgs"){
    stop('only lbfgs is implemented')
  }
  n = length(x)
  opt = try(optim(init_val,
                  fn = vga_binom_obj,
                  gr = vga_binom_obj_grad,
                  x=x,
                  nb=nb,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  gh_points=gh_points,
                  method = 'L-BFGS-B'),silent = TRUE)

  if(class(opt)=='try-error'){
    warning(paste(opt[1],'; returning the initialization values'))
    return(list(m=init_val[1:n],v=exp(init_val[(n+1):(2*n)])))
  }else{
    return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
  }


}

vga_binomial_lbfgs = function(init_val,x,nb,beta,sigma2,method='lbfgs',gh_points){
  if(method!="lbfgs"){
    stop('only lbfgs is implemented')
  }
  n = length(x)
  opt = try(lbfgs::lbfgs(
    call_eval = vga_binom_obj,
                  call_grad = vga_binom_obj_grad,
                  vars=init_val,
                  x=x,
                  nb=nb,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  gh_points=gh_points),silent = TRUE)

  if(class(opt)=='try-error'){
    warning(paste(opt[1],'; returning the initialization values'))
    return(list(m=init_val[1:n],v=exp(init_val[(n+1):(2*n)])))
  }else{
    return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
  }


}

vga_binom_obj = function(theta,x,nb,beta,sigma2,n,gh_points){
  m = theta[1:n]
  lv = theta[(n+1):(2*n)]
  val = - sum(x*m - nb*Elog1pexp(m,exp(lv),gh_points) - (m^2+exp(lv)-2*m*beta)/2/sigma2 + lv/2)
  return(val)
}

vga_binom_obj_grad = function(theta,x,nb,beta,sigma2,n,gh_points){
  m = theta[1:n]
  lv = theta[(n+1):(2*n)]
  dm = - (x - nb*Elog1pexp_dm(m,exp(lv),gh_points) - (m-beta)/sigma2)
  dlv = - (- nb*Elog1pexp_dlv(m,lv,gh_points) - exp(lv)/2/sigma2 + 1/2)
  return(c(dm,dlv))
}

Elog1pexp = function(m,v,gh_points){
  mat = outer(sqrt(2*v),gh_points$x,FUN='*') + m
  return(c(vebpm:::log1pexp(mat)%*%gh_points$w)/sqrt(pi))
}

Elog1pexp_dm = function(m,v,gh_points){
  mat = outer(sqrt(2*v),gh_points$x,FUN='*') + m
  return(c(sigmoid(mat)%*%gh_points$w)/sqrt(pi))
}

Elog1pexp_dlv = function(m,lv,gh_points){
  temp = outer(sqrt(2*exp(lv)),gh_points$x,FUN='*')
  mat = temp + m
  mat2 = temp/2
  return(c((mat2*sigmoid(mat))%*%gh_points$w)/sqrt(pi))
}

# Elog1pexp = function(m,v,gh_points){
#   res=c()
#   for(i in 1:length(m)){
#     res[i] = ghQuad(log1pexp_func,gh_points,m=m[i],v=v[i])/sqrt(pi)
#   }
#   return(res)
# }
#
# log1pexp_func = function(x,m,v){
#   return(vebpm:::log1pexp(sqrt(2*v)*x+m))
# }

logit = function(x){
  log(x/(1-x))
}

sigmoid = function(x){
  1/(1+exp(-x))
}

Esigmoid = function(m,v,gh_points){
  mat = outer(sqrt(2*v),gh_points$x,FUN='*') + m
  return(c(sigmoid(mat)%*%gh_points$w)/sqrt(pi))
}

# Esigmoid = function(m,v,gh_points){
#   res=c()
#   for(i in 1:length(m)){
#     res[i] = ghQuad(sigmoid_func,gh_points,m=m[i],v=v[i])/sqrt(pi)
#   }
#   return(res)
# }
#
# sigmoid_func = function(x,m,v){
#   return(sigmoid(sqrt(2*v)*x+m))
# }


binomial_mean_GG = function(x,nb,beta=NULL,sigma2=NULL,
                            est_prior_mean =TRUE,
                            est_prior_var = TRUE,
                            m_init = NULL,
                            v_init = NULL,
                            maxiter=100,tol=1e-5,n_gh = 10,printevery = 1){

  n = length(x)
  if(is.null(m_init)){
    m = logit(x/nb)
    m[m==-Inf] = logit(0.1)
    m[m==Inf] = logit(0.9)
  }else{
    m = m_init
  }

  if(is.null(v_init)){
    v = rep(1/n,n)
  }else{
    v = v_init
  }

  gh_points = gaussHermiteData(n_gh)
  const = sum(lfactorial(nb)-lfactorial(x)-lfactorial(nb-x))
  obj = -Inf

  for(iter in 1:maxiter){
    if(est_prior_mean){
      beta = mean(m)
    }
    if(est_prior_var){
      sigma2 = mean(m^2+v-2*m*beta+beta^2)
    }


    opt = vga_binomial(c(m,log(v)),x,nb,beta,sigma2,gh_points=gh_points)
    m = opt$m
    v = opt$v

    # ELBO
    obj[iter+1] = sum(x*m-nb*Elog1pexp(m,v,gh_points)) +const - n/2*log(2*pi*sigma2)
    - sum(m^2 + v + beta^2 - 2*m*beta)/2/sigma2 + sum(log(2*pi*v))/2 - n/2

    if(iter%%printevery==0){
      print(paste('At iter', iter, 'elbo=',round(obj[iter+1],3)))
    }
    if((obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      if((obj[iter+1]-obj[iter])<0){
        warning('An iteration decreases ELBO. This is likely due to numerical issues.')
      }
      break
    }
  }

  return(list(posterior = list(mean_logit = m,
                               var_logit = v,
                               mean = Esigmoid(m,v,gh_points)),
              fitted_g = list(beta=beta,sigma2=sigma2),
              elbo=obj[length(obj)],
              obj_trace = obj
              ))

}


binomial_smooth = function(x,nb,sigma2_init=NULL, est_sigma2=TRUE,maxiter=100,tol=1e-4,n_gh=10,filter.number = 1,
                           Eb_init = NULL,
                           family = 'DaubExPhase',ebnm_params=list(mode=0)){
  n = length(x)
  W = (t(GenW(n,filter.number,family)))[-1,]
  obj = -Inf
  sigma2 = sigma2_init
  gh_points = gaussHermiteData(n_gh)
  #Eb = logit(mean(x))
  if(is.null(Eb_init)){
    Eb = rep(0,n)
  }else{
    Eb = Eb_init
  }
  opt = vga_binomial(c(Eb,rep(1/n,n)),x,nb,Eb,sigma2,gh_points=gh_points)
  m = opt$m
  v = opt$v
  for(iter in 1:maxiter){

    qb = smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
    Eb = qb$posterior$mean
    Eb2 = qb$posterior$var + Eb^2

    opt = vga_binomial(c(m,log(v)),x,nb,Eb,sigma2,gh_points=gh_points)
    m = opt$m
    v = opt$v

    if(est_sigma2){
      sigma2 = mean(m^2+v+Eb2-2*m*Eb)
    }

  }
  return(list(Eb=Eb,m=m,sigma2=sigma2))
}

# n = 512
# nb = rep(1,n)
# mu = c(rep(-3,n/2),rep(3,n/2))
# p = sigmoid(mu)
# x = rbinom(n,nb,p)
# plot(x,col='grey80')
# lines(sigmoid(mu))
# lines(sigmoid(fit$Eb),col=4)



