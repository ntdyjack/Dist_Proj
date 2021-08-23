#qrsh -l mem_free=20G -l h_vmem=20G -l h_fsize=20G -pe local 8
library(rhdf5)
library(HDF5Array)
library(mbkmeans)
library(TENxBrainData)
library(BiocParallel)

k <- 15
batch <- 500
sce <- TENxBrainData()
set.seed(123)
keep <- which(DelayedArray::rowSums(counts(sce)) != 0 )
res <- mbkmeans(x=counts(sce[keep,]),clusters=k, batch_size = batch, num_init=1, max_iters=100,verbose=T,compute_labels=F)

one_centroid <- function(x, data) {
    DelayedArray::colSums((x - data)^2)
}

all_labels <- function(data, centroids) {
    ss <- apply(centroids, 1, one_centroid, data = data)
    apply(ss, 1, which.min)
}

labels <- unlist(blockApply(counts(sce[keep,]), all_labels, centroids = res$centroids,
  BPPARAM = BiocParallel::MulticoreParam(10),grid = colAutoGrid(counts(sce[keep,]))))



writeLines(as.character(unname(labels)),"/users/ndyjack/Dist_Proj/tables/test_datasets/1M_neurons_labels.txt")
output_file <- writeHDF5Array(counts(sce[keep,]),
  filepath="/users/ndyjack/Dist_Proj/tables/test_datasets/1M_neurons_full_expr.h5",
  name="10xBD_Full",H5type="H5T_STD_I32LE")

#library(TENxBrainData)
#library(mbkmeans)
#library(BiocParallel)
#sce <- TENxBrainData()
#k <- 15
#batch <- 500

#res_1kcell <- mbkmeans(x=counts(sce[,1:1e3]),clusters=k, batch_size = batch, num_init=1, max_iters=100,verbose=F,compute_labels=T)
#length(res_1kcell$Clusters)

#res_5kcell <- mbkmeans(x=counts(sce[,1:5e3]),clusters=k, batch_size = batch, num_init=1, max_iters=100,verbose=F,compute_labels=T)
#length(res_5kcell$Clusters)

#labels_5k <- blockApply(counts(sce[,1:5e3]), all_labels, centroids = res_5kcell$centroids, BPPARAM = BiocParallel::MulticoreParam(8))
