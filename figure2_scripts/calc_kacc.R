#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
p_load(Matrix,RcppCNPy,reticulate,pbapply)
np <-import("numpy")

#args <- c("/users/ndyjack/Dist_Proj/pickles/","/users/ndyjack/Dist_Proj/tables/test_datasets/","jurkat293t_rnaseq_2cl_hvg500","0","30")
k <- as.integer(args[5])

#cell cluster identity vector
#get rid of hvg string
#cluster_lab <- as.integer(readLines(paste0(args[2],sub("_[h|f].*","",args[3]),"_labs.txt")))
cluster_lab <- as.integer(as.factor(readLines(paste0(sub('pickles','labels',args[1]),sub('hvg1k.mtx','labs.txt',args[3])))))
n <- length(cluster_lab)
ind_mat <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))

#read in cell distance vector
#dist_vec <- readRDS(paste0(args[1],args[3],".",args[4],".py"))
dist_mat <- np$load(paste0(args[1],args[3],".",args[4],".npy"))
#dist_mat <- matrix(0,nrow=n,ncol=n)
#dist_mat[upper.tri(dist_mat)] <- dist_vec
#dist_mat <- t(dist_mat)
#dist_mat[upper.tri(dist_mat)] <- dist_vec
#rm(dist_vec) #memory cleanup

cat("\nestimating kacc...\n")
kacc_vec <- sapply(1:n, function(i) {
  knn_tmp <- order(dist_mat[,i],decreasing=F)[2:(k+1)]
  ret_tmp <- sum(ind_mat[knn_tmp,i])/k
  return(ret_tmp)
})
kacc_fin <- mean(kacc_vec)

#knn_mat <- apply(dist_mat,2,function(y) order(y,decreasing=F)[2:(k+1)])
#for each cell, take find its K closest neighbors
#for each of these knn, compute % intra-cluster
#take mean
#n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
#Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
#dist_inter <- sort(dist_vec[ind_inter])
#dist_intra  <- sort(dist_vec[ind_intra])
#rm(dist_vec) #memory cleanup
#orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
#orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
#gplus_samp <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nz_samp
cat("\nfinished kacc...\n")

output <- paste0(args[3],".",args[4]," ",kacc_fin)
outfile <- paste0(sub(x=args[2],pattern="test_datasets","kacc_new"),args[3],".",args[4],".k",args[5],".txt")
cat("\n writing output to ",outfile,"\n")
write(output,file=outfile)
