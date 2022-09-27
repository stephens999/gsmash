# set.seed(100)
# n = 100
# p = n
# X = matrix(0,nrow=n,ncol=n)
# for(i in 1:n){
#   X[i:n,i] = 1
# }
#
# btrue = rep(0,n)
# btrue[50] = 3
# y = rpois(n,exp(X %*% btrue + rnorm(n)))
# plot(y)
# lines(exp(X %*% btrue))
#
#
# fit.smash = smashr::smash.poiss(y)
#
# lines(fit.smash)
#


# Simulate a data set.
set.seed(1)
n          <- 1000
p          <- 300
X          <- matrix(runif(n*p),n,p)
beta       <- double(p)
beta[1:3] <- c(1,-2,3)
sigma2 = 0.25
y          <- drop(X %*% beta + sqrt(sigma2)*rnorm(n))
plot(y)
y = rpois(n,exp(y))
plot(y)

library(ebnm)
temp = ebnm(beta,rep(1e-10,p),prior_family ='normal_scale_mixture')
temp$fitted_g
#sa2 = c(1e-6, 1e-5, 1e-4, 4e-4, 1e-3, 4e-3, 1e-2, 2e-2, 4e-2, 8e-2, 0.16, 0.32, 0.64, 1, 2, 4, 8, 16)
sa2 = c(1e-3,5.29)
w = c(0.99,0.01)

fit = poisson_nugget_glm_GMGM(X,y,
                             beta.init = beta,
                             sa2=sa2,w=w,sigma2=sigma2,control = list(update.sigma2=FALSE,update.pi=FALSE),printevery = 1)





fit = poisson_nugget_glm_GMG(X,y,
                             beta.init = beta,
                             sa2=sa2,w=w,sigma2=sigma2,control = list(update.sigma2=FALSE,update.pi=FALSE))
plot(fit$obj,type='l')
plot(fit$m)
fit$m[1:5]
fit$Phi[1:5,]
plot(fit$e)


G_m(update_m(X,y,v,e,d,Phi,Sa2,X2,m),X,y,v,e,d,Phi,Sa2,X2)
G_v(update_v(X,y,m,e,d,Phi,Sa2,X2,v),X,X2,m,e,d,Phi,Sa2)
G_e(update_e(X,y,m,v,d,sigma2,X2,e),X,y,m,v,d,sigma2,X2)
G_d(update_d(X,y,m,v,e,sigma2,X2,d),X,X2,m,v,e,sigma2)
