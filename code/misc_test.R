library(fastTopics)
library(ebpmf)
library(Matrix)
data(pbmc_facs)

## original data fit
counts = pbmc_facs$counts
counts_filtered = counts[,colSums(counts!=0)>10]
dim(counts_filtered)

### study the init sigma2 on model fitting.

fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-8),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e8.rds')
rm(fit0)
gc()
fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-6),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e6.rds')
rm(fit0)
gc()
fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-4),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e4.rds')
rm(fit0)
gc()
fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-2),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e2.rds')
rm(fit0)
gc()
fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-1),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e1.rds')
rm(fit0)
gc()
fit0 = ebpmf_log(counts_filtered,
                 general_control = list(maxiter=100,save_init_val=TRUE),
                 init_control = list(single_gene_ebpm = FALSE,n_cores=5,init_tol=1e-1,conv_type='sigma2abs'),
                 flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                 sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(fit0,'output/pbmc_fasttopics/nonnegLF_pe_inittol1e1_sigma2abs.rds')
rm(fit0)
gc()

# pbmc3k_sparse = ebpmf_log(counts_filtered,flash_control = list())
# saveRDS(pbmc3k_sparse,'output/pbmc_fasttopics/pbmc3k_sparse.rds')
pbmc3k_nonnegL = ebpmf_log(counts_filtered,
                           flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_normal),loadings_sign = 1),
                           general_control = list(save_init_val=TRUE),
                           init_control = list(n_cores=5))
saveRDS(pbmc3k_nonnegL,'output/pbmc_fasttopics/pbmc_nonnegL_pe.rds')
# pbmc3k_nonnegLF = ebpmf_log(counts_filtered,
#                             flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
#                             general_control = list(maxiter=30,save_init_val=TRUE,save_latent_M=TRUE),
#                             vga_control = list(n_cores=5),
#                             sigma2_control=list(return_sigma2_trace=TRUE,cap_var_mean_ratio=0))

pbmc3k_nonnegLF = ebpmf_log(counts_filtered,
                            general_control = list(maxiter=30),
                            init_control = list(M_init = pbmc3k_nonnegL$init_val$M_init,sigma2_init = pbmc3k_nonnegL$init_val$sigma2_init),
                            flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),factors_sign=1,loadings_sign = 1),
                            sigma2_control=list(return_sigma2_trace=TRUE))
saveRDS(pbmc3k_nonnegLF,'output/pbmc_fasttopics/pbmc_nonnegLF_pe.rds')

# pbmc3k_nonnegL = ebpmf_log(counts_filtered,flash_control = list(ebnm.fn=c(ebnm::ebnm_unimodal_nonnegative,ebnm::ebnm_point_normal),loadings_sign = 1))
# saveRDS(pbmc3k_nonnegL,'output/pbmc_fasttopics/pbmc3k_nonnegL_un.rds')
# pbmc3k_nonnegLF = ebpmf_log(counts_filtered,flash_control = list(ebnm.fn=c(ebnm::ebnm_unimodal_nonnegative,ebnm::ebnm_unimodal_nonnegative),factors_sign=1,loadings_sign = 1))
# saveRDS(pbmc3k_nonnegLF,'output/pbmc_fasttopics/pbmc3k_nonnegLF_un.rds')

# library(stm)
# library(vebpm)
# l0 = log(rowSums(counts_filtered))
# f0 = log(colSums(counts_filtered)/sum(exp(l0)))
# n = nrow(counts_filtered)
# p = ncol(counts_filtered)
# p0 = 500
# M = matrix(nrow=n,ncol=p0)
# sigma2_init = rep(0,p0)
# for(j in 1:p0){
#   if(T){
#     if(j%%50==0){
#       cat(paste(j,'...'))
#     }
#   }
#   fit = suppressWarnings(ebpm_normal(counts_filtered[,j],g_init = list(mean=l0+f0[j],var=NULL),
#                                      fix_g = c(T,F),tol=1e-5))
#   M[,j] = fit$posterior$mean_log
#   sigma2_init[j] = fit$fitted_g$var
# }
