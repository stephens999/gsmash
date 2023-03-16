
source('code/poisson_STM/plot_factors.R')
source('code/poisson_STM/plot_factors_general.R')
source('code/poisson_STM/structure_plot.R')
source('code/poisson_STM/get_loadings_order.R')
source('code/poisson_STM/order_topic.R')
###################### simulation plot #############################

res = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc3k_simulation/simu_pbmc_fasttopics_large_true_var_large_init_var.rds')
loadings_order = get_loadings_order(res$sim_data$Loading,res$sim_data$Factor,
                                    grouping = pbmc_facs$samples$subpop,n_samples = 5000)
plot0=structure_plot_general(res$sim_data$Loading,res$sim_data$Factor,
                             show_legend = F,
                             grouping =pbmc_facs$samples$subpop,title = 'True',print_plot = F,
                             loadings_order = loadings_order)
i = 3
order_nmf = order_topic(res$sim_data$Loading[,-c(1,2)],res$output[[i]]$fitted_model$nmf$W[,-1])
plot1 = structure_plot_general((res$output[[i]]$fitted_model$nmf$W[,-1])[,order_nmf],
                               (t(res$output[[i]]$fitted_model$nmf$H[-1,]))[,order_nmf],
                               grouping=pbmc_facs$samples$subpop,
                               title='log transformation+NMF(squared loss)',
                               print_plot = F,
                               show_legend = F,
                               loadings_order=loadings_order,
                               remove_l0f0 = F)
order_ebpmf = order_topic(res$sim_data$Loading,res$output[[i]]$fitted_model$ebpmf$fit_flash$L.pm)
plot2 = structure_plot_general(res$output[[i]]$fitted_model$ebpmf$fit_flash$L.pm[,order_ebpmf],
                               res$output[[i]]$fitted_model$ebpmf$fit_flash$F.pm[,order_ebpmf],
                               grouping =pbmc_facs$samples$subpop,
                               title='EBPMF',
                               print_plot = F,
                               show_legend = F,
                               loadings_order=loadings_order)
order_flash = order_topic(res$sim_data$Loading[,-c(1,2)],res$output[[i]]$fitted_model$flash$L.pm[,-1])
plot3 = structure_plot_general((res$output[[i]]$fitted_model$flash$L.pm[,-1])[,order_flash],
                               (res$output[[i]]$fitted_model$flash$L.pm[,-1])[,order_flash],
                               grouping =pbmc_facs$samples$subpop,
                               title='log transformation+EBNMF',
                               remove_l0f0 = F,
                               print_plot = F,
                               show_legend = F,
                               loadings_order=loadings_order)
grid.arrange(plot0,plot2,plot3,plot1,nrow=4)






# true pbmc data


library(fastTopics)
library(Matrix)
data("pbmc_facs")
counts = pbmc_facs$counts
counts = counts[,colSums(counts!=0)>10]
sum(counts==0)/prod(dim(counts))
dim(counts)


# rework cell types

cell_names = pbmc_facs$samples$celltype
cell_names = as.character(cell_names)
cell_names[cell_names=='CD8+ Cytotoxic T'] = 'CD8+ T cells'
cell_names[cell_names%in%c('CD4+ T Helper2','CD4+/CD45RO+ Memory','CD8+/CD45RA+ Naive Cytotoxic','CD4+/CD45RA+/CD25- Naive T','CD4+/CD25 T Reg')] = 'Other T cells'
table(cell_names)
cell_names = as.factor(cell_names)
## plot ebpmf, flash, nmf fit on the same plot
fit = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc_fasttopics/init_tol_effect/nonnegLF_pe_inittol1e2_iter60.rds')
plot_ebpmf = structure_plot_general(fit$fit_flash$L.pm,
                                    fit$fit_flash$F.pm,
                                    #pbmc_facs$samples$celltype,
                                    #pbmc_facs$samples$subpop,
                                    n_samples = 5000,
                                    cell_names,
                                    show_legend = T,
                                    title = 'EBPMF')

library(randomcoloR)
n <- 20
set.seed(12345)
palette <- distinctColorPalette(n)
fit_flash = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc_fasttopics/othermethods/flash_nonnegLF.rds')
kset = order(fit_flash$pve,decreasing = T)
plot_flash = structure_plot_general(fit_flash$L.pm[,kset[2:13]],
                                    fit_flash$F.pm[,kset[2:13]],
                                    #pbmc_facs$samples$celltype,
                                    #pbmc_facs$samples$subpop,
                                    n_samples = 5000,
                                    cell_names,
                                    #colors = palette,
                                    remove_l0f0 = F,
                                    show_legend = F,
                                    title='log transformation + EBNMF')

## plot nmf

# library(NNLM)
# s = rowSums(counts)
# y_tilde = log(1+counts/s*median(s)/0.5)
# fit_nmf = NNLM::nnmf(as.matrix(y_tilde),k=20,method='lee',loss='mse')
fit_nmf = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc_fasttopics/othermethods/nmf_on_log.rds')
temp = ebpmf:::poisson_to_multinom(t(fit_nmf$H),fit_nmf$W)
plot_nmf = structure_plot_general(temp$L[,-1],temp$FF[,1],
                                  #pbmc_facs$samples$celltype,
                                  #pbmc_facs$samples$subpop,
                                  cell_names,
                                  #colors = palette,
                                  n_samples=5000,
                                  remove_l0f0 = F,
                                  show_legend = F,
                                  topic_model = T,
                                  title = 'log transformation + NMF(squared loss)')

grid.arrange(plot_ebpmf,plot_flash,plot_nmf,nrow=3)

## run flash on M
# get M
n = nrow(counts)
p = ncol(counts)
M = ebpmf:::vga_pois_solver_mat_newton(tcrossprod(fit$fit_flash$L.pm,fit$fit_flash$F.pm),
                                       as.matrix(counts),
                                       1,
                                       tcrossprod(fit$fit_flash$L.pm,fit$fit_flash$F.pm),
                                       matrix(fit$sigma2,nrow=n,ncol=p,byrow = T),
                                       maxiter = 1000,
                                       return_V = F)

# remove l0,f0
M0 = M - tcrossprod(fit$fit_flash$L.pm[,c(1,2)],fit$fit_flash$F.pm[,c(1,2)])
library(flashier)
flash_on_M = try(flash.init(M0,var.type = 2) %>%
                              flash.set.verbose(1)%>%
                              flash.add.greedy(Kmax = 20,
                                               ebnm.fn = c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),
                                               init.fn =function(f) init.fn.default(f, dim.signs = c(1, 1))) %>%
                              flash.backfit() %>%
                              flash.nullcheck())

structure_plot_general(flash_on_M$L.pm,
                       flash_on_M$F.pm,
                       pbmc_facs$samples$celltype,
                       colors = palette,remove_l0f0 = F,
                       show_legend = F,
                       title='EBNMF on pseudo data from EBPMF')

### run UMAP

library(dplyr)
library(Seurat)
library(patchwork)

# Initialize the Seurat object with the raw (non-normalized data).
rownames(counts) = pbmc_facs$samples$celltype
pbmc <- CreateSeuratObject(counts = t(counts), project = "pbmc3k", min.cells = 3, min.features = 100,names.field=1)
pbmc <- NormalizeData(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
Idents(pbmc) = pbmc_facs$samples$subpop
DimPlot(pbmc, reduction = "umap")

######3

hist(M - tcrossprod(fit$fit_flash$L.pm[,1],fit$fit_flash$F.pm[,1]),breaks = 100)


########
# [1] "GNLY"   "CCL5"   "NKG7"   "GZMB"   "FGFBP2" "GZMH"   "CLIC3"
# [8] "GZMA"   "PRF1"   "CST7"   "TYROBP" "KLRB1"  "FCGR3A" "CCL4"
# [15] "FCER1G" "KLRF1"  "SPON2"  "PTGDS"  "KLRD1"  "KLRC1"
f1 = sort(fit$fit_flash$F.pm[,3],decreasing = T)[1:2000]
cols = rep('grey80',length(f1))
cols[1:4] = 'red'
plot(f1,xlab = '',
     pch=20,ylab='',col=cols,
     axes = F)

axis(2,at=c(0,1,2,3),labels = c(0,1,2,3))
text(110,f1[1],label='GNLY',cex=0.8)
text(110,f1[2],label='CCL5',cex=0.8)
text(110,f1[3],label='NKG7',cex=0.8)
text(110,f1[4],label='GZMB',cex=0.8)


########
gene_names = pbmc_facs$genes$symbol[which(colSums(pbmc_facs$counts!=0)>10)]
cols = rep('grey80',length(gene_names))
cols[which(gene_names=='GNLY')] = 'red'
cols[which(gene_names=='CCL5')] = 'red'
cols[which(gene_names=='NKG7')] = 'red'
cols[which(gene_names=='GZMB')] = 'red'
cols[which(gene_names=='FGFBP2')] = 'red'
plot(fit$fit_flash$F.pm[,3],xlab = '',
     pch=20,ylab='',col=cols,
     axes = F)
axis(2,at=c(0,1,2,3),labels = c(0,1,2,3))
text(which(gene_names=='GNLY')+600,f1[1],label='GNLY',cex=0.8)
text(which(gene_names=='CCL5')+600,f1[2],label='CCL5',cex=0.8)
text(which(gene_names=='NKG7')+600,f1[3],label='NKG7',cex=0.8)
text(which(gene_names=='GZMB')+600,f1[4],label='GZMB',cex=0.8)
text(which(gene_names=='FGFBP2')+700,f1[5],label='FGFBP2',cex=0.8)
