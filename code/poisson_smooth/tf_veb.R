
library(smashr)
library(genlasso)
library(glmgen)
library(vebpm)
tf_veb = function(x,k,sigma2,s2,
                  maxiter=1000,
                  tol=1e-5,
                  mu_init=NULL,
                  est_s2=FALSE,
                  ebnm_params=list(mode=0,prior_family='normal_scale_mixture')){

  n = length(x)
  obj = rep(0,maxiter+1)
  obj[1] = -Inf
  # init mu
  if(is.null(mu_init)){
    Emu = ti.thresh(x,sigma2)
  }else{
    Emu = mu_init
  }

  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }

  D = getDtf(n,k-1)
  DtD = crossprod(D)
  p = nrow(D)
  d = tfMultiply(Emu,k=k)
  for(iter in 1:maxiter){
    # update b

    a = ebnm(d,sqrt(s2*sigma2),
             mode=ebnm_params$mode,
             prior_family=ebnm_params$prior_family,
             scale = ebnm_params$scale,
             g_init = ebnm_params$g_init,
             fix_g = ebnm_params$fix_g,
             output = ebnm_params$output,
             optmethod = ebnm_params$optmethod)
    Eb = a$posterior$mean
    Vb = a$posterior$sd^2
    dKL = a$log_likelihood + p*(log(2*pi*sigma2*s2)/2)+sum((d^2-2*d*Eb+Eb^2+Vb)/s2/sigma2/2)


    # update mu
    Vmu = chol2inv(chol(diag(n)+DtD/s2))
    Emu = drop(Vmu%*%(x+crossprod(D,Eb)/s2))

    d = tfMultiply(Emu,k=k)
    #
    if(est_s2){
      s2 = (sum(d^2)+sum(diag(DtD%*%Vmu))-2*sum(Eb*d)+sum(Eb^2+Vb))/p
    }

    # ELBO
    obj[iter+1] = tf_veb_obj(x,s2,sigma2,Emu,Vmu,Eb,Vb,dKL,D,DtD)
    print(paste("Done iter",iter,"obj =",obj[iter+1]))
    if(abs(obj[iter+1]-obj[iter])<tol){
      obj = obj[1:(iter+1)]
      break
    }

  }
  return(list(Emu=Emu,Eb=Eb,Vmu=Vmu,Vb=Vb,obj_value=obj,s2=s2,ebnm_res = a))

}

tf_veb_obj = function(x,s2,sigma2,Emu,Vmu,Eb,Vb,dKL,D,DtD){
  d = tfMultiply(Emu,k=k)
  n = length(x)
  return(-(sum(d^2)+sum(diag(Vmu))-2*sum(x*Emu))/2/sigma2-((sum(d^2)+sum(diag(DtD%*%Vmu))-2*sum(Eb*d)+sum(Eb^2+Vb)))/s2/sigma2/2-n*log(s2)/2+determinant(Vmu,logarithm =T)$modulus/2 + dKL)
}





# /**
#   * @brief Multiplies a vector by D transpose, without having
# * to explictly construct or use the matrix D.
# *
#   * @param x                    locations of the responses
# * @param n                    number of observations
# * @param k                    order of the trendfilter
# * @param a                    the input vector to multiply
# * @param b                    allocated space for the output
# * @return void
# * @see tf_dtxtil
# */
#   void tf_dtx(double *x, int n, int k, double *a, double *b)
# {
#   int i;
#   int j;
#   double fact;
#
#   for(i=0; i < n-k; i++) b[i] = a[i];
#
#   if( k < 1 || k >= n )
#     return;
#
#   for(i=k; i > 0; --i)
#   {
#
#     /* b[0:n-i] = D' * b[0:n-i-1] for 1 <= i < n */
#     b[n-i] = b[n-i-1];
#     for(j=n-i-1; j > 0; --j)
#     {
#       b[j] = b[j-1] - b[j];
#     }
#     b[0] = -b[0];
#
#     if( i != 1 )
#     {
#       for(j=0; j <= n-i; ++j)
#       {
#         b[j] = b[j] / ( x[j+i-1] - x[j] );
#       }
#     }
#   }
#
#   fact = glmgen_factorial(k-1);
#   for(i=0; i < n; ++i)
#   {
#     b[i] *= fact;
#   }
# }


