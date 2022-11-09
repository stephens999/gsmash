get_summary_mean = function(out,rm_method = NULL){
  mse_all = c()
  for(i in 1:length(out$output)){
    mse_all = rbind(mse_all,out$output[[i]]$MSE_mean)
  }
  method_names = colnames(mse_all)
  if(!is.null(rm_method)){
    rm_idx = match(rm_method,method_names)
    mse_all = mse_all[,-rm_idx]
  }
  mse_mean = apply(mse_all,2,mean,na.rm=T)
  mse_sd = apply(mse_all,2,sd,na.rm=T)
  order_idx = order(mse_mean)
  res = cbind(mse_mean[order_idx],mse_sd[order_idx])
  colnames(res) = c('mean','sd')
  print(knitr::kable(round(res,3)))
  # datax = reshape2::melt(mse_all,varnames =c('simu','methods'),value.name = 'mse')
  # print(datax %>%
  #         ggplot(aes(x=methods, y=mse)) +
  #         geom_boxplot()+
  #         coord_flip())

  #relative to mle mse
  mse_mle = apply((out$sim_data$X-out$sim_data$Mean)^2,1,mean)
  mse_relative = mse_all/mse_mle
  datax = reshape2::melt(mse_relative,varnames =c('simu','methods'),value.name = 'mse_relative_to_mle')
  print(datax %>%
          ggplot(aes(x=methods, y=mse_relative_to_mle)) +
          geom_boxplot()+
          coord_flip()+
          geom_hline(yintercept = 1, linetype="dashed",
                     color = "red"))
}

get_summary_mean_log = function(out,rm_method = NULL){
  mse_all = c()
  for(i in 1:length(out$output)){
    mse_all = rbind(mse_all,out$output[[i]]$MSE_log_mean)
  }
  method_names = colnames(mse_all)
  if(!is.null(rm_method)){
    rm_idx = match(rm_method,method_names)
    mse_all = mse_all[,-rm_idx]
  }
  mse_mean = apply(mse_all,2,mean,na.rm=T)
  mse_sd = apply(mse_all,2,sd,na.rm=T)
  order_idx = order(mse_mean)
  res = cbind(mse_mean[order_idx],mse_sd[order_idx])
  colnames(res) = c('mean','sd')
  print(knitr::kable(round(res,3)))
  # datax = reshape2::melt(mse_all,varnames =c('simu','methods'),value.name = 'mse_log_mean')
  # print(datax %>%
  #         ggplot(aes(x=methods, y=mse_log_mean)) +
  #         geom_boxplot()+
  #         coord_flip())

  #relative to mle mse
  mse_mle = apply((log(out$sim_data$X+1)-out$sim_data$log_Mean)^2,1,mean)
  mse_relative = mse_all/mse_mle
  datax = reshape2::melt(mse_relative,varnames =c('simu','methods'),value.name = 'mse_relative_to_mle')
  print(datax %>%
          ggplot(aes(x=methods, y=mse_relative_to_mle)) +
          geom_boxplot()+
          coord_flip()+
          geom_hline(yintercept = 1, linetype="dashed",
                     color = "red"))
}
