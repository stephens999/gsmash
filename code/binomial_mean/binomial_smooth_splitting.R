#'@title binomial smoothing via splitting
#'
#'

binom_smooth_split = function(x,nb,
                              n_gh = 10,
                              m_init = NULL,
                              sigma2_init = NULL,
                              maxiter=100,
                              tol=1e-5,
                              filter.number = 1,
                              family = 'DaubExPhase',
                              verbose=TRUE,
                              printevery = 1,
                              ebnm_params=list(mode=0),
                              convergence_criteria = 'objabs',
                              W=NULL,
                              est_sigma2 = TRUE,
                              warmstart = TRUE,
                              make_power_of_2 = 'reflect'){
  t_start = Sys.time()
  n = length(x)
  if(length(nb)==1){
    nb = rep(nb,n)
  }
  m = m_init
  sigma2=sigma2_init
  v = rep(sigma2/2,n)

  sigma2_trace = c(sigma2)
  gh_points = gaussHermiteData(n_gh)

  qb = list(fitted_g = NULL)

  const = sum(lfactorial(nb)-lfactorial(x)-lfactorial(nb-x))
  obj = -Inf

  for(iter in 1:maxiter){
    if(warmstart){
      qb = suppressWarnings(smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=list(g_init=qb$fitted_g),W=W))
    }else{
      qb = smash_dwt(m,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
    }

    Eb = qb$posterior$mean
    Eb2 = qb$posterior$var + Eb^2


    # print(paste('elbo after q(b) is',binom_smooth_split_obj(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,qb$dKL)))
    # browser()

    opt = vga_binomial(c(m,log(v)),x,nb,Eb,sigma2,gh_points=gh_points)
    m = opt$m
    v = opt$v

    # print(paste('elbo after q(mu) is',binom_smooth_split_obj(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,qb$dKL)))

    if(est_sigma2){
      sigma2_new = mean(m^2+v+Eb2-2*m*Eb)
      sigma2_trace = c(sigma2_trace,sigma2_new)
      sigma2 = sigma2_new
    }

    # print(paste('elbo after sigma2 is',binom_smooth_split_obj(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,qb$dKL)))

    obj[iter+1] = binom_smooth_split_obj(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,qb$dKL)
    if(verbose){
      if(iter%%printevery==0){
        print(paste("Done iter",iter,"obj =",round(obj[iter+1],3)))
      }
    }

    if(abs(obj[iter+1]-obj[iter])/n <tol){
      break
    }
  }
  return(list(posterior=list(mean=sigmoid(m),
                             mean_logit = m,
                             mean_smooth = sigmoid(Eb),
                             mean_logit_smooth=Eb),
              fitted_g = list(sigma2=sigma2,sigma2_trace=sigma2_trace),
              elbo=obj[length(obj)],
              elbo_trace = obj,
              H = (qb$dKL + sum(log(2*pi*v)/2-log(2*pi*sigma2)/2-(m^2+v-2*m*Eb+Eb2)/2/sigma2)),
              log_likelihood = obj[length(obj)],
              run_time = difftime(Sys.time(),t_start,units='secs')))
}

# binom_smooth_split_obj = function(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,H){
#   n = length(x)
#   sum(x*m-nb*Elog1pexp(m,v,gh_points)) +const - n/2*log(2*pi*sigma2)
#   - sum(m^2 + v + Eb2 - 2*m*Eb)/2/sigma2 + H + sum(log(2*pi*v))/2 - n/2
# }

binom_smooth_split_obj = function(x,nb,sigma2,m,v,Eb,Eb2,gh_points,const,H){
  n = length(x)
  term1 = sum(x*m-nb*Elog1pexp(m,v,gh_points)-(m^2 + v + Eb2 - 2*m*Eb)/2/sigma2 + log(v)/2)
  term2 = const - n/2*log(2*pi*sigma2)+ H  - n/2
  return(term1 + term2)
}


