

library(flashr)
library(parallel)
library(vebpm)
splitting_PMF = function(Y,S,sigma2=NULL,est_sigma2 = TRUE,
                         Kmax=10,var_type='by_row',
                         M_init = NULL,maxiter=100,tol=0.01,
                         n_cores = 10,verbose_flash=FALSE,
                         printevery=10){

  n = nrow(Y)
  p = ncol(Y)

  if(is.null(S)){
    S = 1
  }
  if(length(S)==1){
    S = matrix(S,nrow=n,ncol=p)
  }
  if(is.null(M_init)){
    M_init = log((Y+1)/S)
  }
  M = M_init
  if(is.null(sigma2)){
    temp = flash(M,Kmax=Kmax,var_type = var_type,verbose = FALSE)
    sigma2 = 1/temp$fit$tau
  }
  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=n,ncol=p)
  }

  const = sum(Y*log(S)) - sum(lfactorial(Y))

  V = matrix(1/n,nrow=n,ncol=p)
  obj = -Inf

  for(iter in 1:maxiter){

    # solve flash
    fit_flash = flash(flash_set_data(M,S=sqrt(sigma2)), Kmax=Kmax,var_type = 'zero',verbose = verbose_flash)
    KL_LF = sum(unlist(fit_flash$fit$KL_f)) + sum(unlist(fit_flash$fit$KL_l))

    #print(paste('After flash,elbo is',calc_split_PMF_obj(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))


    #browser()
    # solve vgap
    res = mclapply(1:n,function(i){
      opt = vga_pois_solver(M[i,],Y[i,],S[i,],fit_flash$fitted_values[i,],sigma2[i,])
      return(opt)
    },mc.cores = n_cores)
    M = do.call(rbind,lapply(res,function(z){z$m}))
    V = do.call(rbind,lapply(res,function(z){z$v}))

    #print(paste('After vga,elbo is',calc_split_PMF_obj(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))

    # update sigma2
    if(est_sigma2){
      R2 = flashr:::flash_get_R2(flash_set_data(M),fit_flash$fit)
      if(var_type=='constant'){
        sigma2 = mean(R2+V)
        sigma2 = matrix(sigma2,nrow=n,ncol=p)
      }else if(var_type=='by_row'){
        sigma2 = rowMeans(V+R2)
        sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = F)
      }else if(var_type=='by_col'){
        sigma2 = colMeans(V+R2)
        sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = T)
      }else{
        stop('Non-supported var type')
      }
    }

    #print(paste('After sigma2,elbo is',calc_split_PMF_obj(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))

    # check convergence
    obj[iter + 1] = calc_split_PMF_obj(Y,S,sigma2,M,V,fit_flash,KL_LF,const)
    if((obj[iter+1] - obj[iter]) < tol){
      if((obj[iter+1] - obj[iter])<0){
        warning('An iteration decreases ELBO')
      }
      break
    }

    if(iter%%printevery==0){
      print(paste('At iter',iter, 'ELBO=',obj[iter+1]))
    }

  }

  return(list(M=M,V=V,fit_flash=fit_flash,elbo=obj[length(obj)],eblo_trace=obj,sigma2 = sigma2))
}



calc_split_PMF_obj = function(Y,S,sigma2,M,V,fit_flash,KL_LF,const){
  sum(Y*M - S*exp(M+V/2) - log(2*pi*sigma2)/2 - (V+flashr:::flash_get_R2(flash_set_data(Y=M),fit_flash$fit))/2/sigma2 + log(2*pi*V)/2 + 0.5 ) + const+ KL_LF
}













