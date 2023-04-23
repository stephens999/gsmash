#####
library(readr)
R.utils::gunzip("/project2/mstephens/dongyue/poisson_mf/droplet/GSE103354_Trachea_droplet_UMIcounts.txt.gz")
suppressMessages(
  counts <- read_delim('/project2/mstephens/dongyue/poisson_mf/droplet/GSE103354_Trachea_droplet_UMIcounts.txt',
                       delim = "\t",progress = FALSE,
                       col_types = cols(.default = col_double(),
                                        gene = col_character())))


dat = readRDS("/project2/mstephens/dongyue/poisson_mf/droplet_66k/pulseseq.rds")
dim(dat)
dat[1:6,1:6]
devs <- scry::devianceFeatureSelection(dat)
dev_ranked_genes <- rownames(dat)[order(devs, decreasing = TRUE)]
topdev <- head(dev_ranked_genes, 3000)
'Krt5'%in%topdev
'Ascl1'%in%topdev
'Ascl2'%in%topdev
'Foxi1'%in%topdev
'Foxi1'%in%dev_ranked_genes
'Foxj1'%in%topdev
'Cdhr3'%in%topdev
'Cdhr3'%in%dev_ranked_genes
'Atp6v0d2'%in%topdev
'Chga'%in%topdev
which(dev_ranked_genes=='Chga')

gene_idx = rownames(dat)%in%topdev
summary(rowSums(dat[gene_idx,]))
gene_names = genes$symbol[gene_idx]


library(Seurat)
droplet <- CreateSeuratObject(counts = dat, project = "droplet", min.cells = 3, min.features = 200)
droplet <- NormalizeData(droplet)
droplet <- FindVariableFeatures(droplet, selection.method = "vst", nfeatures = 3000)
topdev_seurat = VariableFeatures(droplet)
'Chga'%in%topdev_seurat
'Foxi1'%in%topdev_seurat
'Chga'%in%topdev_seurat
'Ascl1'%in%topdev_seurat
'Ascl2'%in%topdev_seurat
'Ascl3'%in%topdev_seurat
'Rgs13'%in%topdev_seurat
'Cdhr3'%in%topdev_seurat
'Foxj1'%in%topdev_seurat
'Krt5'%in%topdev_seurat
'Krt17'%in%topdev_seurat
'Cftr'%in%topdev_seurat
'Atp6v0d2'%in%topdev_seurat


selected_genes = union(topdev,topdev_seurat)

dat = dat[rownames(dat)%in%c(selected_genes),]
# filter out genes have expression in less than 10 cells
dat = dat[rowSums(dat!=0)>10,]
dim(dat)
saveRDS(t(dat),file='/project2/mstephens/dongyue/poisson_mf/droplet_66k/pulse_seq_5kgene.rds')

cell.type <- sapply(strsplit(colnames(dat), "_"), `[[`, 5)
time.point <- sapply(strsplit(colnames(dat), "_"), `[[`, 2)

cell.type <- ifelse(
  cell.type == "Basal", paste0("Basal (", time.point, ")"), cell.type
)

cell.type <- factor(cell.type, levels = c(
  "Basal (Tp0)",
  "Basal (Tp30)",
  "Basal (Tp60)",
  "Club",
  "Club (hillock-associated)",
  "Ciliated",
  "Goblet.progenitor",
  "Goblet.1",
  "Goblet.2",
  "Tuft.progenitor",
  "Tuft.1",
  "Tuft.2",
  "Neuroendocrine",
  "Ionocyte",
  "Proliferating"
))

#################
######liger######
liger.dat <- rliger::createLiger(list(dat = dat))
liger.dat <- rliger::normalize(liger.dat)
liger.dat <- rliger::selectGenes(liger.dat)


