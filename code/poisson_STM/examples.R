

set.seed(12345)
N = 100
p = 50
K = 2
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1
Ftrue[21:40,2] = 1
Ltrue = matrix(rnorm(N*K), ncol=K)

l0 = runif(N,1,2)
f0 = runif(p,1,2)
S0 = tcrossprod(l0,f0)

Lambda = S0*exp(tcrossprod(Ltrue,Ftrue))

Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)

fit = splitting_PMF(Y,S0,Kmax = 100)
plot(fit$obj,type='l')
fit$fit_flash$nfactors

plot(fit$fit_flash$fitted_values,tcrossprod(Ltrue,Ftrue))

plot(fit$fit_flash$ldf$f[,1])
plot(fit$fit_flash$ldf$f[,2])
plot(fit$fit_flash$ldf$f[,3])
plot(fit$fit_flash$ldf$f[,4])


fit$sigma2[,1]

X = log((Y+1)/S0)


#X = matrix(rnorm(N*p,mean=tcrossprod(Ltrue,Ftrue)),nrow=N,ncol=p)
fit1 = flash(X,var_type = 'constant')
fit1$fit$tau
fit1$objective

fit2 = flash(flash_set_data(Y=X,S=sqrt(1/fit1$fit$tau)),var_type = 'zero')
fit2$objective
