#options("width"=160)
#library(Seurat)
library(PoiClaClu)
library(scran)
library(DropletUtils)

#load in null data
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19") 
counts <- as.matrix(counts(sce))
counts <- counts[-which(rowSums(counts)==0),]

#compute poisson distance
dist.pois <- PoissonDistance(t(counts))$dd

#normalize with scran
sce <- SingleCellExperiment(list(counts=counts))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
counts.norm <- logcounts(sce)
#compute euclidean distance on the scran-normalized counts
dist.euc <- dist(t(counts.norm))

save.image("/users/ndyjack/Dist_Proj/rdata/cd14_noclec7a.rda")
saveRDS(dist.euc,"/users/ndyjack/Dist_Proj/rdata/cd14_euc.rda/")
saveRDS(dist.euc,"/users/ndyjack/Dist_Proj/rdata/cd14_poi.rda/")
