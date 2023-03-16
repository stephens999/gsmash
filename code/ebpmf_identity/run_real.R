datax = read.csv('/project2/mstephens/gtex-stm/Counts/TPM3.Counts.csv.gz')
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
fit_stm = ebpmf_identity(datax,K=7,init = list(L_init = fit_fasttopics$L,F_init = fit_fasttopics$F))
plot(fit_stm$res$qf$Esmooth_f[,1],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,2],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,3],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,4],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,5],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,6],col='grey80',pch='.',cex=2)
plot(fit_stm$res$qf$Esmooth_f[,7],col='grey80',pch='.',cex=2)

par(mfrow=c(2,1))
for(i in 1:7){
  plot(fit_fasttopics$F[,i],col='grey80',pch='.',cex=2,type='l')
  plot(fit_stm$res$qf$Esmooth_f[,i],col='grey80',pch='.',cex=2,type='l')
}


