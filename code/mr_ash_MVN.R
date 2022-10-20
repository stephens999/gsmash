
library(mr.ash)
library(Rfast)
library(Matrix)

# it's hard to init this method. Especially when sa2 has very small elements.
# For example when min(sa2) = 1e-10, the method breaks entirely. But when min(sa2) = 1e-4, then it works even initialize from NULL.
# so let's use init_mr_ash()
mr_ash_MVN = function(X,y,
                      fit0 = init_mr_ash(X,y),
                      control = list(),
                      printevery = 10,
                      min_sa2 = 1e-4){

  n = length(y)
  p = dim(X)[2]
  # if(is.null(sa2)){
  #   sa2 = (2^((1:20) / 20) - 1)^2
  # }

  control = modifyList(mr_ash_MVN_control_default(),control,keep.null = TRUE)

  # if(is.null(beta.init)){
  #   m = as.vector(double(p))
  # }else{
  #   m = beta.init
  # }
  #
  #
  # if(is.null(w)){
  #   w = rep(1,K)/K
  # }
  m = fit0$b
  sigma2=  fit0$resid.sd^2
  sa2 = fit0$prior$sd^2
  sa2[1] = min_sa2
  K = length(sa2)
  w = fit0$prior$weights
  W = matrix(w,nrow=p,ncol=K,byrow=TRUE)

  # if(is.null(sigma2)){
  #   sigma2 = var(drop(y-X%*%m))
  # }

  #Phi = matrix(rep(w,each = p),nrow = p)
  #V = matrix(0,nrow=p,ncol=p)
  Sa2 = matrix(sa2,nrow=p,ncol=K,byrow=TRUE)
  Phi = log(W) - log(Sa2*sigma2)/2-matrix(m^2,nrow=p,ncol=K,byrow=FALSE)/2/Sa2/sigma2
  Phi = Phi - c(apply(Phi,1,max))
  Phi = exp(Phi)
  Phi = Phi/c(rowSums(Phi))
  a = c(rowSums(Phi/Sa2))
  XtX = Crossprod(X,X)
  Xty = crossprod(X,y)
  y2_sum = sum(y^2)

  obj = c(-Inf)

  for(iter in 1:control$max.iter){

    # update m, V

    U = chol2inv(chol(XtX+diag(a)))
    m = U%*%Xty
    V = sigma2*U

    # update Phi

    Phi = log(W) - log(Sa2*sigma2)/2-matrix(m^2+diag(V),nrow=p,ncol=K,byrow=FALSE)/2/Sa2/sigma2
    Phi = Phi - c(apply(Phi,1,max))
    Phi = exp(Phi)
    Phi = Phi/c(rowSums(Phi))
    Phi = pmax(Phi,sqrt(.Machine$double.eps))

    a = c(rowSums(Phi/Sa2))



    # update w,sigma2
    if(control$update.prior.weights){
      w = colMeans(Phi)
      w = pmax(w,sqrt(.Machine$double.eps))
      W = matrix(w,nrow=p,ncol=K,byrow=TRUE)
    }

    if(control$update.resid.sd){
      sigma2 = c(y2_sum + sum(m*((XtX+diag(a))%*%m)) + sum(diag((XtX+diag(a))%*%V)) - 2*sum(Xty*m))/(n+p)
    }


    # calc obj
    obj[iter+1] = calc_obj_mar_ash_MVN(X,y,m,V,Phi,a,W,sigma2,Sa2,XtX,Xty,y2_sum)
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
  return(list(posteriorMean=m,posteriorVar=V,w=w,sigma2=sigma2,obj=obj,Phi=Phi))

}

mr_ash_MVN_control_default = function(){
  list(min.iter = 1, max.iter = 1000, convtol = 1e-05, update.prior.weights = TRUE,
       update.resid.sd = TRUE, update.order = NULL, epstol = 1e-12)
}

calc_obj_mar_ash_MVN = function(X,y,m,V,Phi,a,W,sigma2,Sa2,XtX,Xty,y2_sum){
  n = length(y)
  p = dim(X)[2]
  K = dim(W)[2]

  val = -n/2*log(sigma2) - c(y2_sum+sum(m*(XtX+diag(a))%*%m) + sum(diag((XtX+diag(a))%*%V)) - 2*sum(Xty*m))/2/sigma2 + as.numeric(determinant(V,log=TRUE)$modulus)/2 + sum(Phi*(log(W)-log(sigma2*Sa2)/2-log(Phi)))

  return(drop(val))
}

simulate_regression_data  = function (n, p, s, sigma = 1, pve = 0.5, center_X = TRUE, standardize_X = FALSE,
          intercept = 0, ncov = 0)
{

  n <- floor(n)
  p <- floor(p)
  s <- floor(s)
  ncov <- floor(ncov)
  if (s > p)
    stop("Input s should be no greater than p")
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = center_X, scale = standardize_X)
  b.idx <- 1:s
  b <- rep(0, p)
  if (s > 0) {
    b.values <- rnorm(s)
    b[b.idx] <- b.values
    if (pve < 1) {
      sb <- sqrt(pve/(1 - pve)/var(drop(X %*% b)))
      b <- sb * sigma * b
    }
  }
  if (ncov > 0) {
    Z <- matrix(rnorm(n * ncov), n, ncov)
    u <- rnorm(ncov)
  }
  else {
    Z <- as.numeric(NA)
    u <- as.numeric(NA)
  }
  y <- intercept + X %*% b + sigma * rnorm(n)
  if (ncov > 0)
    y <- y + Z %*% u
  y <- drop(y)
  rnames <- paste0("s", 1:n)
  cnames_X <- paste0("v", 1:p)
  names(y) <- rnames
  rownames(X) <- rnames
  colnames(X) <- cnames_X
  names(b) <- cnames_X
  if (ncov > 0) {
    cnames_Z <- paste0("c", 1:ncov)
    names(u) <- cnames_Z
    colnames(Z) <- cnames_Z
    rownames(Z) <- rnames
  }
  return(list(y = y, X = X, Z = Z, u = u, b = b, intercept = intercept,
              sigma = sigma))
}

tf_basis = function(n,k){
  # H = np.zeros((n, n))
  # npowerk = np.power(n, k)
  # seq = np.arange(1, n+1).reshape(n, 1)
  # H[:, :k + 1] = np.power(np.tile(seq, k+1), np.arange(k+1)) / np.power(n, np.arange(k+1))
  # for j in range(k+1, n):
  #     for i in range(n):
  #         if i > j - 1:
  #             Hij = 1.0
  #             for l in range(1, k+1):
  #                 Hij *= (i - j + k - l + 1)
  #             H[i, j] = Hij #/ np.power(n, k)
  H = matrix(0,nrow=n,ncol=n)

}

