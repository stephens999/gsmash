

library(smashrgen)
library(parallel)
source("~/Rpackages/gsmash/code/poisson_smooth/simu_study_poisson_smooth.R")
snr_list = c(1,3)
count_size_list = c(5,10,100)
smooth_func_list = c('sblocks','cblocks','heavi','angles','bursts','spike')
n=1024
n_simu = 30
output_path = 'output/poisson_smooth_simulation/'
for(smooth_func in smooth_func_list){
  for(snr in snr_list){
    for(count_size in count_size_list){
      print(paste('Running:',smooth_func, 'snr:',snr,'count-size:',count_size))
      datax = sim_data_smooth(n_simu=n_simu,n=n,snr=snr,
                              count_size = count_size,
                              smooth_func = smooth_func)
      out = simu_study_poisson_smooth(datax,n_cores=10)
      saveRDS(out,file=paste(output_path,smooth_func,n_simu,'_n_',n,'_snr_',snr,'_count_size_',count_size,'.rds',sep=''))
    }
  }
}


######## test

datax = sim_data_smooth(n_simu=2,n=128,snr=1,
                        count_size = 10,
                        smooth_func = 'sblocks')
out = simu_study_poisson_smooth(datax,n_cores=5)




#
# out <- readRDS("~/Documents/projPhD/gsmash/output/poisson_smooth_simulation/bumps30_n_512_snr_3_count_size_100.rds")
#
# plot(out$sim_data$X[1,],col='grey80')
# lines(exp(out$sim_data$b),col='grey60')
# lines(out$output[[1]]$fitted_model$smash$posterior$mean_smooth,col='2')
# lines(out$output[[1]]$fitted_model$split_ndwt$posterior$mean_smooth,col='2')
# lines(out$output[[1]]$fitted_model$smash_two_step$posterior$mean_smooth,col='2')
#
#
# temp = smashrgen:::pois_smooth_split(out$sim_data$X[1,],Eb_init = out$sim_data$X[1,],wave_trans = 'dwt',sigma2_init = 0.1,est_sigma2 = F)
# plot(out$sim_data$X[1,],col='grey80')
# lines(exp(out$sim_data$b),col='grey60')
# lines(temp$posterior$mean_smooth,col=4)
# temp$fitted_g
