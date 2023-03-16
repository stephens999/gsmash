# library(fastTopics)
# library(flashier)
# library(glmpca)
# library(Matrix)
# load('/project2/mstephens/pcarbo/git/single-cell-topics/data/droplet.RData')
# K = 20
# glmpca_droplet_poi = glmpca(counts,L = K, fam='poi',minibatch = 'memoized')
# saveRDS(glmpca_droplet_poi,file='/project2/mstephens/dongyue/poisson_mf/droplet/glmpca_droplet_poi.rds')
#
# glmpca_droplet_nb = glmpca(counts,L = K, fam='nb',minibatch = 'memoized')
# saveRDS(glmpca_droplet_nb,file='/project2/mstephens/dongyue/poisson_mf/droplet/glmpca_droplet_nb.rds')
#
# s = rowSums(counts)
# y_tilde = log(1+counts/s*median(s)/0.5)
# y_tilde=Matrix(y_tilde,sparse=T)
#
# flash_droplet = flash(y_tilde,var.type = 2,greedy.Kmax = K,backfit = T)
# saveRDS(flash_droplet,file='/project2/mstephens/dongyue/poisson_mf/droplet/flash_droplet.rds')
#
# flash_droplet_nonnegL = flash(y_tilde,var.type = 2,greedy.Kmax = K,backfit = T,ebnm.fn = c(ebnm::ebnm_unimodal_nonnegative,ebnm::ebnm_point_normal))
# saveRDS(flash_droplet_nonnegL,file='/project2/mstephens/dongyue/poisson_mf/droplet/flash_droplet_nonnegL.rds')
#
# flash_droplet_nonnegLF = flash(y_tilde,var.type = 2,greedy.Kmax = K,backfit = T,ebnm.fn = c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential))
# saveRDS(flash_droplet_nonnegLF,file='/project2/mstephens/dongyue/poisson_mf/droplet/flash_droplet_nonnegLF.rds')
#



library(Matrix)
res = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc3k_simulation/simu_pbmc_fasttopics_large_true_var_large_init_var.rds')
for(i in c(1,2,3,4,5)){
  S = rowSums(res$sim_data$Y[[i]])
  Y_normed = log(1+median(S)*res$sim_data$Y[[i]]/S/0.5)
  res$output[[i]]$fitted_model$flash = try(flash.init(Y_normed,var.type = 2) %>%
                                             flash.set.verbose(1)%>%
                                             flash.add.greedy(Kmax = 9,ebnm.fn = c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),
                                                              init.fn =function(f) init.fn.default(f, dim.signs = c(1, 1))) %>%
                                             flash.backfit() %>%
                                             flash.nullcheck())
  res$output[[i]]$fitted_model$nmf = NNLM::nnmf(as.matrix(Y_normed),k=9,method='lee',loss='mse')
}

saveRDS(res,'/project2/mstephens/dongyue/poisson_mf/pbmc3k_simulation/simu_pbmc_fasttopics_large_true_var_large_init_var.rds')
source('code/poisson_STM/get_loadings_order.R')
loadings_order = get_loadings_order(res$sim_data$Loading,res$sim_data$Factor,
                                    grouping = pbmc_facs$samples$subpop,n_samples = 5000)
plot0=structure_plot_general(res$sim_data$Loading,res$sim_data$Factor,
                             grouping =pbmc_facs$samples$subpop,title = 'True',print_plot = F,
                             loadings_order = loadings_order)
i = 1
plot1 = structure_plot_general(res$output[[i]]$fitted_model$nmf$W[,-1],
                             t(res$output[[i]]$fitted_model$nmf$H[-1,]),
                             grouping=pbmc_facs$samples$subpop,
                             title='log transformation+NMF(squared loss)',
                             print_plot = T,
                             loadings_order=loadings_order,
                             remove_l0f0 = F)
  plot2 = structure_plot_general(res$output[[i]]$fitted_model$ebpmf$fit_flash$L.pm,
                                 res$output[[i]]$fitted_model$ebpmf$fit_flash$F.pm,
                                 grouping =pbmc_facs$samples$subpop,
                                 title='ebpmf',
                                 print_plot = F,
                                 loadings_order=loadings_order)
  plot3 = structure_plot_general(res$output[[i]]$fitted_model$flash$L.pm[,-1],
                                 res$output[[i]]$fitted_model$flash$L.pm[,-1],
                                 grouping =pbmc_facs$samples$subpop,
                                 title='log transformation+flash',
                                 remove_l0f0 = F,
                                 print_plot = T,
                                 loadings_order=loadings_order)
  grid.arrange(plot0,plot1,plot2, plot3,nrow=4)



