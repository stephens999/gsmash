ridge_obj = function(lambda,y,A,v){
  n = length(y)
  return(-mgcv::dmvn(y,rep(0,n),1/lambda*A+diag(v)))
}

ebs = function(x,y,v,k){
  #tk <- seq(min(x),max(x),length=k) ## knots
  #X<-apply(diag(k),1,function(z) approx(tk,z,x,rule=2)$y)
  X = bs(x,df=k)
  D <- diff(diag(k),d=2)
  S <- crossprod(D)
  A = X%*%S%*%t(X)
  lbd = optimize(ridge_obj,c(0.01,10000),y=y,A=A,v=v)
  llk = lbd$objective
  lbd = lbd$minimum
  b_post = solve(crossprod(X)+lbd*S)%*%t(X)%*%y
  return(list(fitted.value = X%*%b_post,loglik = -llk,lambda=lbd))

}

n = 500
#v = runif(n,0.5,1.5)
v = rep(1,n)
#y = c(rep(0,n/4),rep(3,n/4),rep(0,n/4),rep(2,n/4)) + rnorm(n,0,sqrt(v))
x = seq(0,1,length.out = n)
y = sin(2*pi*x) * 2 + rnorm(n,0,sqrt(v))


fit = ebs(x,y,v,k=50)

plot(x,y,col='grey80',pch=19)
lines(x,fit$fitted.value)
fit$loglik
fit$lambda

k = 50
X = bs(x,df=k)
D <- diff(diag(k),d=2)
S <- crossprod(D)
A = X%*%S%*%t(X)
lls = seq(0.01,100,length.out = 100)
oo = c()
for(i in 1:length(lls)){
  oo[i] = ridge_obj(lls[i],y,A,v)
}
plot(lls,oo,type='l')


fitg = gam(y ~ s(x,bs = "ps",k=k), family=gaussian,method = 'ML',scale=1,weights = 1/v)
lines(x,fitg$fitted.values,col=4)
