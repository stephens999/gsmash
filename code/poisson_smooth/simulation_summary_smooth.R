get_mse_mean = function(out,rm_method=NULL){
  res = c()
  for(i in 1:length(out$output)){
    res = rbind(res,out$output[[i]]$mse_smooth)
  }
  methods = colnames(res)
  if(!is.null(rm_method)){
    idx = match(rm_method,methods)
    res = res[,-idx]
    methods = methods[-idx]
  }
  n_methods = length(methods)

  res_mean = colMeans(res,na.rm = T)
  res_sd = apply(res,2,sd,na.rm=T)

  order_idx = order(res_mean,decreasing = F)
  print(knitr::kable(cbind(res_mean[order_idx],res_sd[order_idx]),col.names = c('mean','sd')))

  #mse_mle = apply((out$sim_data$X-out$sim_data$Mean)^2,1,mean)
  #mse_relative = mse_all/mse_mle
  datax = reshape2::melt(res,varnames =c('simu','methods'),value.name = 'mse')
  print(datax %>%
          ggplot(aes(x=methods, y=mse)) +
          geom_boxplot()+
          coord_flip())

}


get_runtime = function(out,rm_method=NULL){
  res = c()
  for(i in 1:length(out$output)){
    res = rbind(res,out$output[[i]]$mse_smooth)
  }

  res_time = c()
  for(i in 1:length(out$output)){
    res_time = rbind(res_time,out$output[[i]]$run_times)
  }

  methods = colnames(res)
  if(!is.null(rm_method)){
    idx = match(rm_method,methods)
    res = res[,-idx]
    res_time = res_time[,-idx]
    methods = methods[-idx]
  }
  n_methods = length(methods)

  res_mean = colMeans(res,na.rm=T)
  res_time_mean = colMeans(res_time,na.rm=T)

  plot(log2(res_time_mean),res_mean,xlab='run time(log2)',ylab='mse',pch = 1:n_methods,col=1:n_methods, main='mean estimation')
  legend('topright',methods,pch=1:n_methods,col=1:n_methods)
  #abline(h=1,lty=2,col='grey80')

}



plot_all_curves = function(out,method='smash_two_step_hetero'){
  #x_min = apply(out$sim_data$X,2,min)
  #x_max = apply(out$sim_data$X,2,max)
  x_median = apply(out$sim_data$X,2,median)
  plot(x_median,col='grey80',xlab='',ylab='x(median of all runs)',main=method)
  #lines(x_min,type='p',col='grey80')
  lines(exp(out$sim_data$b))
  for(i in 1:length(out$output)){
    if(class(out$output[[i]]$fitted_model[[method]])!='try-error'){
      lines(out$output[[i]]$fitted_model[[method]]$posterior$mean_smooth,col='grey80')
    }

  }


  res = c()
  for(i in 1:length(out$output)){
    res = rbind(res,out$output[[i]]$mse_smooth)
  }
  methods = colnames(res)
  res0 = res[,match(method,methods)]
  lines(out$output[[which.min(res0)]]$fitted_model[[method]]$posterior$mean_smooth,col=2,lwd=1.5)
  lines(out$output[[which.max(res0)]]$fitted_model[[method]]$posterior$mean_smooth,col=4,lwd=1.5)

  legend('topright',c('True mean','best fit','worst fit','all others'),lty=c(1,1,1,1),col=c(1,2,4,'grey80'))


}












