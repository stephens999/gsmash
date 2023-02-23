source('code/ebpmf_simulation/simu_func.R')

# load real data fit results

fit = readRDS("~/Documents/myproj/gsmash/output/pbmc_fasttopics/pbmc3k_nonnegLF_pe.rds")
simdata = sim_data_pmf_log(fit$fit_flash$L.pm,fit$fit_flash$F.pm,fit$sigma2,var_type = 'by_col',n_simu=3)
rm(fit)
gc()

res = simu_study_PMF(simdata,n_cores = 3,
                     method_list=c('flash','ebpmf'),
                     ebnm_function = ebnm::ebnm_point_exponential,
                     loadings_sign = 1,
                     factors_sign = 1,
                     Kmax=10,
                     var_type='by_col',
                     maxiter=100,
                     tol=1e-5)
saveRDS(res,file='output/ebpmf_simulation/pbmc_fasttopics_simu.rds')
