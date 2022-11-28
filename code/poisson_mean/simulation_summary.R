get_summary_runtime = function(out,rm_method = NULL,include_log_res = TRUE,log_runtime = TRUE){
  runtime = c()
  for(i in 1:length(out$output)){
    runtime = rbind(runtime,out$output[[i]]$run_times)
  }
  method_names = colnames(runtime)
  if(!is.null(rm_method)){
    rm_idx = match(rm_method,method_names)
    runtime = runtime[,-rm_idx]
  }
  runtime_mean = apply(runtime,2,median,na.rm=T)
  runtime_sd = apply(runtime,2,sd,na.rm=T)
  order_idx = order(runtime_mean)
  res = cbind(runtime_mean[order_idx],runtime_sd[order_idx])
  colnames(res) = c('median','sd')
  print(knitr::kable(round(res,3)))



  # plot runtime vs mse

  mse_all = c()
  for(i in 1:length(out$output)){
    mse_all = rbind(mse_all,out$output[[i]]$MSE_mean)
  }
  method_names = colnames(mse_all)
  if(!is.null(rm_method)){
    rm_idx = match(rm_method,method_names)
    mse_all = mse_all[,-rm_idx]
    method_names = method_names[-rm_idx]
  }
  n_method = length(method_names)
  mse_mle = apply((out$sim_data$X-out$sim_data$Mean)^2,1,mean)
  mse_relative = mse_all/mse_mle
  mse_mean = apply(mse_relative,2,mean,na.rm=T)
  #mse_sd = apply(mse_all,2,sd,na.rm=T)

  if(log_runtime){
    plot(log2(runtime_mean),mse_mean,xlab='run time(log2)',ylab='mse relative to mle',pch = 1:n_method,col=1:n_method, main='mean estimation')
    legend('topright',method_names,pch=1:n_method,col=1:n_method)
    abline(h=1,lty=2,col='grey80')
  }else{
    plot(runtime_mean,mse_mean,xlab='run time',ylab='mse relative to mle',pch = 1:n_method,col=1:n_method, main='mean estimation')
    legend('topright',method_names,pch=1:n_method,col=1:n_method)
    abline(h=1,lty=2,col='grey80')
  }


  if(include_log_res){
    # plot runtime vs mse_log
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
    if(log_runtime){
      plot(log2(runtime_mean),mse_mean,xlab='run time(log2)',ylab='mse',pch = 1:n_method,col=1:n_method, main='log mean estimation')
      legend('topright',method_names,pch=1:n_method,col=1:n_method)
    }else{
      plot(runtime_mean,mse_mean,xlab='run time',ylab='mse',pch = 1:n_method,col=1:n_method, main='log mean estimation')
      legend('topright',method_names,pch=1:n_method,col=1:n_method)
    }

  }


}

get_summary_elbo = function(out,rm_method = c( "nb_lb","nb_pg","ash_pois_log")){
  elbo_list = lapply(out$output,function(x){
    lapply(x$fitted_model,function(z){
      if('elbo'%in%names(z)){
        return(z$elbo)
      }else if('loglik'%in%names(z)){
        return(z$loglik)
      }else if('log_likelihood'%in%names(z)){
        return(z$log_likelihood)
      }else{
        return(NA)
      }
    })
  })
  elbo_list = do.call(rbind,lapply(elbo_list,unlist))

  method_names = colnames(elbo_list)
  if(!is.null(rm_method)){
    rm_idx = match(rm_method,method_names)
    elbo_list = elbo_list[,-rm_idx]
  }

  elbo_mean = apply(elbo_list,2,mean,na.rm=T)
  elbo_sd = apply(elbo_list,2,sd,na.rm=T)
  order_idx = order(elbo_mean,decreasing = T)
  res = cbind(elbo_mean[order_idx],elbo_sd[order_idx])
  colnames(res) = c('mean','sd')
  print(knitr::kable(round(res,3)))


}

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

get_summary_mean_log = function(out,rm_method = NULL,rm_x0 = TRUE){
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

  datax = reshape2::melt(mse_all,varnames =c('simu','methods'),value.name = 'mse_log_mean')
  print(datax %>%
          ggplot(aes(x=methods, y=mse_log_mean)) +
          geom_boxplot()+
          coord_flip())

  # #relative to mle mse
  # if(rm_x0){
  #   mse_mle = apply((log(out$sim_data$X+1)-out$sim_data$log_Mean)^2,1,mean)
  #   mse_relative = mse_all/mse_mle
  #   datax = reshape2::melt(mse_relative,varnames =c('simu','methods'),value.name = 'mse_relative_to_mle')
  #   print(datax %>%
  #           ggplot(aes(x=methods, y=mse_relative_to_mle)) +
  #           geom_boxplot()+
  #           coord_flip()+
  #           geom_hline(yintercept = 1, linetype="dashed",
  #                      color = "red"))
  # }else{
  #   mse_mle = apply((log(out$sim_data$X+1)-out$sim_data$log_Mean)^2,1,mean)
  #   mse_relative = mse_all/mse_mle
  #   datax = reshape2::melt(mse_relative,varnames =c('simu','methods'),value.name = 'mse_relative_to_mle')
  #   print(datax %>%
  #           ggplot(aes(x=methods, y=mse_relative_to_mle)) +
  #           geom_boxplot()+
  #           coord_flip()+
  #           geom_hline(yintercept = 1, linetype="dashed",
  #                      color = "red"))
  # }

}
