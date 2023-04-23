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



###########
chr22 = read.table('/project2/compbio/todongyue/gnomad/chr22.expobs.SNV.1kb.bed')
data <- fread('/project2/compbio/todongyue/gnomad/chr22.expobs.SNV.1kb.bed', fill = T)
plot(data$SNV_obs/data$SNV_exp,col='grey80')

data <- data[complete.cases(data$SNV_exp),]
data$SNV_obs[data$SNV_obs ==0] <- 0.5

data$randeff <- data$SNV_obs/data$SNV_exp
data$log_randeff <- log(data$randeff)
data$mean_lgre <- smash.gaus(data$log_randeff)


data$mean_homo = smash.gaus(data$log_randeff,sigma=smashr:::sd_estimate_gasser_etal(data$log_randeff))


plot(data$SNV_obs/data$SNV_exp,col='grey80')
lines(exp(data$mean_homo))
