idx <- seq(1,1000,by=1)
files <- paste0("/fastscratch/myscratch/ndyjack/hashes/tmp_10xBD.",idx,".csv")
hashes <- lapply(files, function(x) read.csv(x,header=F))

zero_ds <- sapply(hashes, function(x) any(rowSums(x)==0))

hashes_combined  <- do.call(rbind,hashes)
cell_labels <- readLines("/users/ndyjack/Dist_Proj/tables/test_datasets/1M_neurons_labels.txt")
zero_hashes <- rowSums(hashes_combined)==0


library(DelayedArray)
setRealizationBackend("HDF5Array")
dat <- HDF5Array("/users/ndyjack/Dist_Proj/tables/test_datasets/1M_neurons_full_expr.h5",name = "10xBD_Full")
lib_sizes <- DelayedArray::colSums(dat)
#sizes_nzh <- DelayedArray::colSums(dat[,which(zero_hashes)])
#sizes_zh <- DelayedArray::colSums(dat[,which(!zero_hashes)])



find_lsh_knn = function(x,dat,k,l){
  #i <- 1
  matches <- lapply(1:l, function(j) {
    which(dat[,j]==x[j])
  })
  kNN_tab <- sort(table(unlist(matches)),decreasing=T)
  kNN <- names(kNN_tab)[1:k] 
 return(kNN)
}

sample_cells <- sample(1:nrow(hashes_combined), size=100)
l_values <- c(10,20,50,100)

kNN_accuracy_by_l <- sapply(l_values, function(l) {
  acc_tmp <- sapply(sample_cells, function(i) {
    #print(i)
    nn_tmp <- find_lsh_knn(x = as.integer(hashes_combined[i,]), dat=hashes_combined[-sample_cells,], k=100,l=l)
    nn_id_tmp <- cell_labels[as.numeric(nn_tmp)]
    id_tmp <- cell_labels[i]
    sum(nn_id_tmp==id_tmp) / length(nn_tmp)
  })
  return(acc_tmp)
})


setRealizationBackend("HDF5Array")
dat <- HDF5Array(paste0(args[1],"1M_neurons_full_expr.h5"),name = "10xBD_Full")
