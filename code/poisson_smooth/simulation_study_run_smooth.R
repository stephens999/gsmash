# snr = 1,3
# count_size = 3, 10, 100
# blocks, bumps, heavi,doppler

library(smashrgen)
snr_list = c(1,3)
count_size_list = c(3,10,100)
smooth_func = "bumps"
n=512
n_simu = 30
output_path = 'output/poisson_smooth_simulation/'
for(snr in snr_list){
  for(count_size in count_size_list){
    print(paste('Running blocks, snr:',snr,'count-size:',count_size))
    datax = sim_data_smooth(n_simu=n_simu,n=n,snr=snr,
                            count_size = count_size,
                            smooth_func = smooth_func)
    out = simu_study_poisson_smooth(datax,n_cores=10)
    saveRDS(out,file=paste(output_path,smooth_func,n_simu,'_n_',n,'_snr_',snr,'_count_size_',count_size,'.rds',sep=''))
  }
}
