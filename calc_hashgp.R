#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")

# args <- c('/users/ndyjack/Dist_Proj/tables/test_datasets/','/users/ndyjack/Dist_Proj/tables/hashes/', 'lsh_testdata_2cl_full','l1')
#cell cluster identity vector
#get rid of hvg string
cluster_lab <- as.integer(readLines(paste0(args[1],sub("_[h|f].*","",args[3]),"_labs.txt")))
Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
Nz <- (Nt*(Nt-1))/2
#cell pair x cluster identity vector
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
tmp <- ind_vec[upper.tri(ind_vec,diag=F)]
ind_intra <- which(tmp)
ind_inter <- which(!tmp)

nhash_vec <-c(seq(10,90,by=20),seq(100,1000,by=300))
#,seq(1000,10000,by=3000))
# c(10,50,100,500,1000,5000,10000) #seq(1000,10000,by=10000) #seq(10,100,by=10))
gplus_vec <- rep(NA,length(nhash_vec))
#read in cell hash matrix
hash_mat <- t(read.csv(paste0(args[2],args[3],".",args[4],".csv"),header = FALSE))

cat("\nestimating g+...\n")
n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
n <- ncol(hash_mat)

for(k in 1:length(nhash_vec)){
  cat("\n",k,"of",length(nhash_vec),"\n")
  #dist_vec <- lapply(1:(n-1), function(i) {
  #  sapply((i+1):(n), function(j) {
  #    #sum(hash_mat[1:nhash_vec[k],i] != hash_mat[1:nhash_vec[k],j]) / nhash_vec[k]
  #    median(abs(hash_mat[1:nhash_vec[k],i] - hash_mat[1:nhash_vec[k],j])) 
  #  })
  #})
  dist_mat <- sapply(1:n, function(i) sapply(1:n, function(j) median(abs(hash_mat[1:nhash_vec[k],i] - hash_mat[1:nhash_vec[k],j]))))
  dist_vec <- dist_mat[upper.tri(dist_mat)]
  dist_inter <- sort(dist_vec[ind_inter])
  dist_intra  <- sort(dist_vec[ind_intra])
  gplus_vec[k] <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nz_samp
}
cat("\nfinished g+...\n")

output_dir <- paste0(sub(x=args[1],pattern="test_datasets","gplus"),args[3],".",args[4],".txt") 
output_file <- cbind(nhash_vec,gplus_vec)
cat("\n writing output to ",output_dir,"\n")
write.csv(x=output_file,file=output_dir,quote=F,row.names=F)
