temp = read.table('/project2/compbio/todongyue/AtoG/mutation.A_to_G.and.mutation.T_to_C.10k')
x = temp$V4
x = as.numeric(x)
x = x[!is.na(x)]
sum(x==0)/length(x)
plot(x[1:100])
fit_tf = trendfilter(x)
fit_tf_cv =cv.trendfilter(fit_tf)
plot(fit_tf,lambda = fit_tf_cv$lambda.min,col='grey80',pch=20)

library(glmgen)
