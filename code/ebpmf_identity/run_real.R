datax = read.csv('/project2/mstephens/gtex-stm/Counts/TPM3.Counts.csv.gz')
rownames(datax) = datax[,1]
datax = datax[,-1]
datax = as.matrix(datax)
datax = datax[,colSums(datax)!=0]
library(Matrix)
library(fastTopics)
X = Matrix(datax,sparse=T)

fit_fasttopics = fit_topic_model(X,k=7)
plot(fit_fasttopics$F[,1],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,2],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,3],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,4],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,5],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,6],col='grey80',pch='.',cex=2)
plot(fit_fasttopics$F[,7],col='grey80',pch='.',cex=2)

structure_plot(fit_fasttopics)

library(ebpmf)
tpm3_fit_fasttopics <- readRDS("~/Rpackages/gsmash/output/tpm3_fit_fasttopics.rds")
fit_stm = ebpmf_identity(datax,K=7,init = list(L_init = tpm3_fit_fasttopics$L,F_init = tpm3_fit_fasttopics$F))
saveRDS(fit_stm,'output/tpm3_fit_stm.rds')
plot(fit_stm$res$qf$Esmooth_f[,1],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,2],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,3],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,4],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,5],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,6],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,7],col='grey80',pch='.',cex=2)

par(mfrow=c(2,1))
for(i in 1:7){
  #plot(fit_fasttopics$F[,i],col='grey50',pch='.',cex=2,type='l')
  plot(tpm3_fit_stm$res$qf$Esmooth_f[,i],col='grey50',pch='.',cex=2,type='l')
  plot(fit_stm$res$qf$Ef_smooth[,i],col='grey50',pch='.',cex=2,type='l')
  # plot(runmed(fit_stm$res$qf$Esmooth_f[,i],k=53),col='grey50',pch='.',cex=2,type='l')
}


sad = c()
k_list = seq(3,333,by=10)
for(k in k_list){
  sad = c(sad,sum(abs(fit_stm$EF[,1] - runmed(fit_stm$EF[,1],k=k))))
}
plot(k_list,sad,type='l')

k = 8
plot(fit_stm$EF[,k],col='grey80',pch='.',cex=3)
lines(runmed(fit_stm$EF[,k],k=43),type='l')
