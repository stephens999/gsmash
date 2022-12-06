library(fastTopics)
library(Matrix)
library(stm)
data(pbmc_facs)
counts <- pbmc_facs$counts
table(pbmc_facs$samples$subpop)
## use only B cell and NK cell
cells = pbmc_facs$samples$subpop%in%c('B cell', 'NK cell')
Y = counts[cells,]
dim(Y)
# filter out genes that has few expressions(3% cells)
genes = (colSums(Y>0) > 0.03*dim(Y)[1])
Y = Y[,genes]
# make sure there is no zero col and row
sum(rowSums(Y)==0)
sum(colSums(Y)==0)
dim(Y)

rm(counts)
rm(genes)
rm(cells)

S = tcrossprod(c(rowSums(Y)),c(colSums(Y)))/sum(Y)
Y = as.matrix(Y)
# run method
fit = splitting_PMF_flashier(Y,S,var_type = 'by_col',Kmax = 30,maxiter = 5000,verbose = TRUE,n_cores = 10)
saveRDS(fit,file='output/poisson_MF_simulation/fit_pbmc_2cells.rds')

# transform Y and fit flashier
a = c(apply(S,2,median))
Y_normed = log((Y/S)*(a/0.5)+1)
Y_normed = Matrix(Y_normed,sparse = T)
fit_flashier = flashier::flash(Y_normed,var.type = 2,greedy.Kmax = 10,backfit = T)
saveRDS(fit_flashier,file='output/poisson_MF_simulation/fit_flashier_pbmc_2cells.rds')

fit_svd = svd(Y_normed)
saveRDS(fit_svd,file='output/poisson_MF_simulation/fit_svd_pbmc_2cells.rds')
