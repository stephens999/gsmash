library(fastTopics)
library(Matrix)
library(stm)
data(pbmc_facs)
counts <- pbmc_facs$counts
table(pbmc_facs$samples$subpop)
## use certain cell types
cells = pbmc_facs$samples$subpop%in%c('B cell', 'NK cell','CD34+')
Y = counts[cells,]
dim(Y)
# filter out genes that has few expressions(3% cells)
genes = (colSums(Y>0) > 0.03*dim(Y)[1])
Y = Y[,genes]
# make sure there is no zero col and row
sum(rowSums(Y)==0)
sum(colSums(Y)==0)
dim(Y)


S = tcrossprod(c(rowSums(Y)),c(colSums(Y)))/sum(Y)
Y = as.matrix(Y)
# run method
fit = splitting_PMF_flashier(Y,S,var_type = 'by_col',Kmax = 30,maxiter = 5000,verbose = TRUE,n_cores = 10,maxiter_backfitting = 1)
saveRDS(fit,file='output/poisson_MF_simulation/fit_pbmc_3cells.rds')

# transform Y and fit flashier
a = c(apply(S,2,median))
Y_normed = log((Y/S)*(a/0.5)+1)
Y_normed = Matrix(Y_normed,sparse = T)
fit_flashier = flashier::flash(Y_normed,var.type = 2,greedy.Kmax = 30,backfit = T,nullcheck = TRUE)
saveRDS(fit_flashier,file='output/poisson_MF_simulation/fit_flashier_pbmc_3cells.rds')

fit_svd = svd(Y_normed)
saveRDS(fit_svd,file='output/poisson_MF_simulation/fit_svd_pbmc_3cells.rds')

########################
counts <- pbmc_facs$counts
S = tcrossprod(c(rowSums(counts)),c(colSums(counts)))/sum(counts)
a = c(apply(S,2,median))
Y_normed = log(t(t((counts/S))*(a/0.5))+1)
rm(S)
Y_normed = Matrix(Y_normed,sparse = T)
fit_flashier = flashier::flash(Y_normed,var.type = 2,greedy.Kmax = 30,backfit = T,nullcheck = TRUE)
saveRDS(fit_flashier,file='output/poisson_MF_simulation/fit_flashier_pbmc.rds')
########################
counts <- pbmc_facs$counts
genes = (colSums(counts>0) > 10)
counts = counts[,genes]
S = tcrossprod(c(rowSums(counts)),c(colSums(counts)))/sum(counts)
a = c(apply(S,2,median))
Y_normed = log(t(t((counts/S))*(a/0.5))+1)
rm(S)
Y_normed = Matrix(Y_normed,sparse = T)
fit_flashier = flashier::flash(Y_normed,var.type = 2,greedy.Kmax = 30,backfit = T,nullcheck = TRUE)
saveRDS(fit_flashier,file='output/poisson_MF_simulation/fit_flashier_pbmc_filter_gene.rds')
