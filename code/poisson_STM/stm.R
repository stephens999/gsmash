#'@title Smoothed Poisson Topic Model
#'@description This function fits Poisson Topic Model with smooth Loading or Factors
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param init_loss loss function of the initialization method, either mkl or mse.
#'@param maxiter maximum iterations
#'@param tol stop criteria
#'@param smooth_method splitting or bmsm
#'@param fix_F if TRUE, F will not be updated.
#'@param bmsm_control control parameters of BMSM, see bmsm_control_default()
#'@param ebpm_method point_gamma or two_gamma
#'@param ebpm_control control parameters of ebpm, see ebpm_control_default()
#'@param splitting_control control parameters of smashgen, see splitting_control_default()
#'@param smooth_f,smooth_l whether to get smooth estimate of loadings or factors.
#'@param nugget whether to assume nugget effects
#'@param return_all whether return all outputs or simplified ones
#'@return EL,EF: posterior of loadings and factors
#'@examples
#'set.seed(123)
#'n = 120
#'p = 256
#'K= 3
#'L = matrix(0, nrow=n, ncol=K)
#'FF = matrix(0, nrow=K, ncol=p)
#'L[1:(n/3),1] = 1
#'L[((n/3)+1):(2*n/3),2] = 1
#'L[((2*n/3)+1):n,3] = 1
#'L = L + matrix(runif(n*K,0,0.5),nrow=n)
#'FF[1,1:(p/3)] = 1+10
#'FF[2,((p/3)+1):(2*p/3)] = 1+10
#'FF[3,((2*p/3)+1):p] = 1+10
#'lambda = L %*% FF
#'X = matrix(rpois(n=length(lambda),lambda),nrow=n)
#'image(X)
#'@import mixsqp
#'@import NNLM
#'@import ebpm
#'@import Matrix
#'@import vebpm
#'@export

stm = function(X,K,
               init = 'scd',
               init_loss = 'mkl',
               maxiter=100,
               tol=1e-3,
               smooth_method = 'splitting',
               fix_F = FALSE,
               bmsm_control=list(),
               splitting_control=list(),
               ebpm_method='normal',
               ebpm_control=list(),
               smooth_f=TRUE,
               smooth_l=FALSE,
               nugget=FALSE,
               printevery=10,
               verbose=TRUE,
               convergence_criteria = 'ELBOabs'){


  #browser()
  n = dim(X)[1]
  p = dim(X)[2]

  res = init_stm(X,K,init,init_loss,smooth_l,smooth_f,nugget)
  EZ = array(dim = c(n,p,K))
  EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

  KL = c()
  KL[1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))

  obj = c()
  obj[1] = -Inf

  # X = Matrix::Matrix(X,sparse = TRUE)
  # X_idx = summary(X)

  if(smooth_f | smooth_l){
    if(!nugget){
      convergence_criteria = 'ELBOabs'
    }
  }

  for(iter in 1:maxiter){
    #b_k_max = 0
    for(k in 1:K){
      # # get row and col sums of EZ_k
      # b_k = res$ql$Elogl[X_idx$i,k]+res$qf$Elogf[X_idx$j,k] - res$a
      # EZ_k = sparseMatrix(i=X_idx$i,j=X_idx$j,x = X_idx$x*exp(b_k)/res$b,dims = c(n,p))
      l_seq = rowSums(EZ[,,k])
      l_scale = sum(res$qf$Ef[,k])

      # Update L
      if(smooth_l){
        lk_hat = update_smooth(l_seq, l_scale, nugget,bmsm_control,splitting_control)
        res$ql$El[,k] = lk_hat$E
        res$ql$Elogl[,k] = lk_hat$Elog
        res$Hl[k] = lk_hat$H
        res$nugget_l[k] = lk_hat$nugget
        res$Esmooth_l = lk_hat$Esmooth
        #res$gl[[k]] = lk_hat$pi_weights

      }else{
        lk_hat = update_nsmooth(l_seq,l_scale,ebpm_control,ebpm_method)
        res$ql$El[,k] = lk_hat$E
        res$ql$Elogl[,k] = lk_hat$Elog
        res$Hl[k] = lk_hat$H
        res$gl[[k]] = lk_hat$fitted_g
      }


      # Update F

      if(!fix_F){
        f_seq = colSums(EZ[,,k])
        f_scale = sum(res$ql$El[,k])

        if(smooth_f){
          fk_hat = update_smooth(f_seq, f_scale,nugget,bmsm_control,splitting_control)
          res$qf$Ef[,k] = fk_hat$E
          res$qf$Elogf[,k] = fk_hat$Elog
          res$Hf[k] = fk_hat$H
          #loglikR = loglikR + fk_hat$loglik
          res$nugget_f[k] = fk_hat$nugget
          res$Esmooth_f = fk_hat$Esmooth
          #res$gf[[k]] = fk_hat$pi_weight

        }else{
          fk_hat = update_nsmooth(f_seq,f_scale,ebpm_control,ebpm_method)
          res$qf$Ef[,k] = fk_hat$E
          res$qf$Elogf[,k] = fk_hat$Elog
          res$Hf[k] = fk_hat$H
          #loglikR = loglikR + fk_hat$log_likelihood
          res$gf[[k]] = fk_hat$fitted_g
        }
      }

      # b_k_new = res$ql$Elogl[X_idx$i,k] + res$qf$Elogf[X_idx$j, k] - res$a
      # res$b = res$b - exp(b_k) + exp(b_k_new)
      # b_k_max = pmax(b_k_new, b_k_max)


    }
    # Update Z
    EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

    # res$b = res$b/exp(b_k_max)
    # res$a = b_k_max + res$a
    if(convergence_criteria == 'mKLabs'){
      KL[iter+1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, mKL: %f',iter,KL[iter+1]))
        }
      }

      if(abs(KL[iter+1]-KL[iter])<=tol){
        break
      }

    }

    if(convergence_criteria=='ELBOabs'){
      obj[iter+1] = calc_stm_obj(X,K,res)
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, ELBO: %f',iter,obj[iter+1]))
        }
      }
      if((obj[iter+1]-obj[iter])<tol){
        break
      }
    }


  }
  if(iter==maxiter){
    warning('Reached maximum iterations')
  }

  #lambda_hat = tcrossprod(res$ql$El,res$qf$Ef)
  #lambda_init = L_init%*%F_init

  # loglik = sum(dpois(X,lambda_hat,log = TRUE))

  if(smooth_l&nugget){
    if(smooth_f&nugget){
      ldf = poisson2multinom(res$qf$Esmooth_f,res$ql$Esmooth_l)
    }else{
      ldf = poisson2multinom(res$qf$Ef,res$ql$Esmooth_l)
    }
  }else{
    if(smooth_f&nugget){
      ldf = poisson2multinom(res$qf$Esmooth_f,res$ql$El)
    }else{
      ldf = poisson2multinom(res$qf$Ef,res$ql$El)
    }
  }

  fit = list(res = res,EL = ldf$L,EF = ldf$FF,d=ldf$s,obj=obj,KL=KL)
  return(fit)
  # if(return_all){
  #   return(list(ql=ql_hat,qf=qf_hat,gf=gf_hat,gl=gl_hat,KL=KL,Lambda_hat=lambda_hat,
  #               init = inited,
  #               input = list(X=X,K=K),nugget=list(nugget_l=nugget_l,nugget_f=nugget_f)))
  # }else{
  #   return(list(ql=ql_hat$El,qf=qf_hat$Ef,nugget=list(nugget_l=nugget_l,nugget_f=nugget_f),KL=KL))
  # }

}

calc_stm_obj = function(X,K,res){
  n = nrow(X)
  p = ncol(X)
  val = 0
  qz = calc_qz(X,K,res$ql,res$qf)
  for(k in 1:K){
    val = val + qz[,,k]*(matrix(res$ql$Elogl[,k],nrow=n,ncol=p,byrow=F)+matrix(res$qf$Elogf[,k],nrow=n,ncol=p,byrow=T)-log(qz[,,k]))
  }
  E1 = sum(X*val) - sum(tcrossprod(res$ql$El,res$qf$Ef))

  return(E1+sum(res$Hl)+sum(res$Hf))
  #return(val-sum(tcrossprod(res$ql$El,res$qf$Ef)) +sum(res$Hl)+sum(res$Hf) )
}

#'@title initialize the stm model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@param init_loss mkl or mse
#'@export
init_stm = function(X,K,init,init_loss,smooth_l,smooth_f,nugget){

  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

    if(is.null(L_init)){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method='lee',
                              loss='mse',show.warning = F,
                              init = list(H=t(F_init)),
                              verbose = F,max.iter = 50)
      L_init = X_init_fit$W
    }

  }else{

    if(init%in%c('scd','lee')){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = 50)
      L_init = X_init_fit$W
      F_init = t(X_init_fit$H)
    }
    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=p,ncol=K)
      ratio = median(X)/(median(L_init)*median(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
    if(init == 'kmeans'){
      kmeans.init=kmeans(as.matrix(X),K,nstart=5)
      L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
      F_init = t(kmeans.init$centers)
      row.names(F_init)=NULL
    }
  }

  # adjust scale of L and F, mainly for stability.
  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init

  Elogl = log(L_init+1e-10)
  Elogf = log(F_init+1e-10)

  if(smooth_l&nugget){
    ql = list(El = L_init, Elogl = Elogl, Esmooth_l=L_init)
  }else{
    ql = list(El = L_init, Elogl = Elogl, Esmooth_l=NULL)
  }
  if(smooth_f&nugget){
    qf = list(Ef = F_init, Elogf = Elogf, Esmooth_f=F_init)
  }else{
    qf = list(Ef = F_init, Elogf = Elogf, Esmooth_f=NULL)
  }


  # a = 0
  # b = 0
  #
  # X = Matrix(X,sparse = TRUE)
  # d = summary(X)
  #
  # temp = Elogl[d$i,] + Elogf[d$j,]
  # a = rowMax(temp)
  # b = rowSums(exp(temp-a))
  #

  gl = list()
  gf = list()

  return(list(ql=ql,qf=qf,gl=gl,gf=gf,
              #a=a,b=b,
              Hl = rep(0,K),
              Hf = rep(0,K),
              nugget_l = ifelse(nugget,rep(0,K),NULL),
              nugget_f = ifelse(nugget,rep(0,K),NULL)))

}


#'@title update L or F, assume they are smooth
#'@param x sequence of observations
#'@param sf scaling factor
#'@param ...

update_smooth = function(x,s,nugget,bmsm_control=list(),splitting_control=list()){
  if(min(x) < 0){stop ("negative values in x not permitted")}
  if(nugget){
    control0 = splitting_control_default()
    if (any(!is.element(names(splitting_control),names(control0))))
      stop("Argument \"splitting_control\" contains unknown parameter names")
    control1 = modifyList(control0,splitting_control,keep.null = TRUE)

    fit = pois_smooth_split(x,s,Eb_init = control1$Eb_init,
                            sigma2_init = control1$sigma2_init,
                            est_sigma2 = control1$est_sigma2,
                            maxiter = control1$maxiter,
                            tol=control1$tol,
                            filter.number = control1$filter.number,
                            family = control1$family,
                            verbose=control1$verbose,
                            printevery = control1$printevery,
                            ebnm_params=control1$ebnm_params,
                            optim_method=control1$optim_method)

    # fit = smash.gen.poiss(x,s=sf,filter.number=control1$filter.number,
    #                       family=control1$family,nugget=control1$nugget,
    #                       robust=control1$robust,
    #                       robust.q = control1$robust.q,
    #                       transformation = control1$transformation,
    #                       method = control1$method,
    #                       nug.init = control1$nug.init,
    #                       ash.pm = control1$ash.pm,
    #                       eps = control1$eps,
    #                       maxiter = control1$maxiter,
    #                       tol = control1$tol)
    est = fit$Emean
    est_log = fit$Emu
    nugget.est = fit$sigma2
    H = fit$H
    results = list("E" = est,
                   'Elog' = est_log,
                   "H" = H,
                   "nugget" = nugget.est,
                   "Esmooth" = exp(fit$Eb))
    return(results)

  }else{
    fit = BMSM(x,s,bmsm_control)
    est = fit$E
    est_log = fit$Elog
    #pi_weights = fit$pi_weights
    results = list("E" = est,
                   'Elog' = est_log,
                   "H" = NULL,
                   "nugget" = NULL)
    return(results)
  }

  #loglik = fit$loglik





}


update_nsmooth = function(x,s,ebpm_control = list(),ebpm_method){

  control0 = ebpm_control_default()
  if (any(!is.element(names(ebpm_control),names(control0))))
    stop("Argument \"ebpm_control\" contains unknown parameter names")
  control1 = modifyList(control0,ebpm_control,keep.null = TRUE)

  #scale = control1$scale
  #point_mass=control1$point_mass
  #nullweight=control1$nullweight
  #shape= control1$shape
  g_init = control1$g_init
  fix_g = control1$fix_g
  #m = control1$m
  control =  control1$control
  #low = control1$low
  #d = control1$d
  pi0 = control1$pi0

  if(ebpm_method=='point_gamma'){
    out = ebpm::ebpm_point_gamma(x,s,g_init,fix_g,pi0,control)
    H = out$log_likelihood - sum(x*log(s)+x*out$posterior$mean_log-out$posterior$mean*s-lfactorial(x))
    #out$H =H
    results = list("E" = out$posterior$mean,
                   'Elog' = out$posterior$mean_log,
                   "fitted_g" = out$fitted_g,
                   "H" = H)
  }

  if(ebpm_method=='normal'){
    out = vebpm::pois_mean_GG(x,s)
    H = sum(-log(2*pi*out$fitted_g$var)/2-(out$fitted_g$mean^2+out$posterior$posteriorMean_log_mean^2 + out$posterior$posteriorVar_log_mean - 2*out$posterior$posteriorMean_log_mean*out$fitted_g$mean)/2/out$fitted_g$var)
    H = H + sum(log(2*pi*out$posterior$posteriorVar_log_mean)/2)
    results = list("E" = out$posterior$posteriorMean_mean,
                   'Elog' = out$posterior$posteriorMean_log_mean,
                   "fitted_g" = out$fitted_g,
                   "H" = H)
  }



  return(results)



}

#'@title Default parameters of ebpm
#'@export
ebpm_control_default = function(){
  list(pi0 = 'estimate',
       g_init = NULL,
       fix_g = FALSE,
       control =  NULL)
}


#'@title Default parameters of smash gen
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details.
#'@export
splitting_control_default = function(){
  list(Eb_init = NULL,
       sigma2_init = NULL,
       est_sigma2 = TRUE,
       maxiter = 100,
       tol=1e-5,
       filter.number = 1,
       family = 'DaubExPhase',
       verbose=FALSE,
       printevery = 10,
       ebnm_params=list(mode=0),
       optim_method='L-BFGS-B')
}


rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}

Calc_EZ = function(X,K,EZ,ql_hat,qf_hat){
  n = nrow(X)
  p = ncol(X)
  for(k in 1:K){
    EZ[,,k] = outer(ql_hat$Elogl[,k], qf_hat$Elogf[,k], "+")
  }
  EZ = softmax3d(EZ)
  EZ = as.vector(EZ)*as.vector(X)
  dim(EZ) = c(n,p,K)
  EZ
}

calc_qz = function(X,K,ql,qf){
  n = nrow(X)
  p = ncol(X)
  qz = array(dim = c(n,p,K))
  for(k in 1:K){
    qz[,,k] = outer(ql$Elogl[,k], qf$Elogf[,k], "+")
  }
  return(softmax3d(qz))
}

softmax3d=function(x){
  #x = x - array(apply(x,c(1,2),max),dim=dim(x))
  x = exp(x)
  p=as.vector(x)/as.vector(rowSums(x,dims=2))
  p = pmax(p,1e-10)
  dim(p) <- dim(x)
  return(p)
}
