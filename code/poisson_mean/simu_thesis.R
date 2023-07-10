library(vebpm)

#################### log link prior #######################
######################3######################3######################3
######################3######################3######################3

n_cores = 10
n = 1000
n_simu=32
prior_var1 = 2
prior_w0_list = c(0)
prior_mean_list=c(5)
method_list = c('GG','GMGM',
                'GMGM_pointmass',
                'nb_pg',
                'log1exp',
                'split',
                'split_mixture',
                'penalty_compound',
                'penalty_inversion',
                'ash_pois_identity',
                'ebpm_gamma',
                'ebpm_exp_mixture')
for(prior_mean in prior_mean_list){
  for(w0 in prior_w0_list){
    datax = gen_data_log_link(n=n,n_simu=n_simu,w=w0,
                              prior_mean=prior_mean,
                              prior_var0 = 0,
                              prior_var1 = prior_var1,
                              seed=12345)
    out = simu_study_poisson_mean(datax,n_cores = n_cores,method_list=method_list)
    saveRDS(out,file=paste('log_link',n_simu,'_n_',n,'_priormean_',prior_mean,'_priorvar1_',prior_var1,'_priorw_',w0*10,'.rds',sep=''))
  }
}




#################### exponential prior #######################
######################3######################3######################3
######################3######################3######################3


library(vebpm)
n_cores = 10
n = 1000
n_simu=32
exp_rate_list = c(0.1,0.5,1)
method_list = c('GG','GMGM',
                'GMGM_pointmass',
                'nb_pg',
                'log1exp','split','split_mixture',
                'penalty_compound',
                'penalty_inversion',
                'ash_pois_identity',
                'ebpm_gamma',
                'ebpm_exp_mixture')
for(exp_rate in exp_rate_list){
  datax = gen_data_exp(n=n,n_simu=n_simu,
                       prior='exponential',
                       exp_rate = exp_rate,
                       seed=12345)
  out = simu_study_poisson_mean(datax,n_cores = n_cores,method_list=method_list)
  saveRDS(out,file=paste('exp_prior',n_simu,'_n_',n,'_mean_',1/exp_rate,'.rds',sep=''))
}

