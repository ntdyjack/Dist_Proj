#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc)

#args <- c("/users/ndyjack/Dist_Proj/rdata/","/users/ndyjack/Dist_Proj/tables/","jurkat293t_rnaseq_2cl_hvg500","0","gplus_results_July2020.txt")
#cell cluster identity vector
#get rid of hvg string
cluster_lab <- as.integer(readLines(paste0(args[2],sub("_[h|f].*","",args[3]),"_labs.txt")))
Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
Nz <- (Nt*(Nt-1))/2
#cell pair x cluster identity vector
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
ind_vec <- apply(ind_vec,1,function(x) {
  tmp <- rep(0,length(x))
  tmp[which(x)] <- cluster_lab[which(x)]
  return(tmp)
})
ind_vec <- ind_vec[upper.tri(ind_vec)]
ind_intra <- which(ind_vec>0)
ind_inter <- which(ind_vec==0)

#read in cell distance vector
dist_vec <- readRDS(paste0(args[1],args[3],".",args[4],".rda"))$'dist'

cat("\nestimating g+...\n")
n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
dist_inter <- sort(dist_vec[ind_inter])
dist_intra  <- sort(dist_vec[ind_intra])
rm(dist_vec) #memory cleanup
orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
gplus_samp <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nz_samp
cat("\nfinished g+...\n")

output <- paste0(args[3],".",args[4]," ",gplus_samp)
outfile <- paste0(sub(x=args[2],pattern="test_datasets","gplus"),args[3],".",args[4],".txt")
cat("\n writing output to ",outfile,"\n")
write(output,file=outfile)
