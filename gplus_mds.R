#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,heatmap3,beeswarm,doParallel)

#cell cluster identity vector
cluster_lab <- as.integer(readLines(paste0(args[1],args[3],"_labs.txt")))
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
dist_vec <- readRDS(paste0(args[2],args[3],".",args[4],".rda"))[[1]]

#parallel backend
#make this more general at some point, just set manually for now

#coerce to full distance matrix
dist_mat <- matrix(0,nrow=length(cluster_lab),ncol=length(cluster_lab))
dist_mat[lower.tri(dist_mat)] <- dist_vec
dist_mat <- t(dist_mat)
dist_mat[lower.tri(dist_mat)] <- dist_vec
dist_mat <- as.dist(dist_mat)

#try with MDS
mds <- cmdscale(d=dist_mat,k=50)
mds <- as.matrix(dist(mds))
dist_vec <- mds[upper.tri(mds)]
#sm <- foreach(k=ind_intra,.combine=c) %dopar% { sum(dist_vec[k] > dist_vec[ind_inter]) }
sm <- mclapply(1:length(ind_intra), function(i) {
  sum(dist_vec[ind_intra[i]] > dist_vec[ind_inter])
}, mc.cores=30,mc.preschedule=T)
sm <- sum(unlist(sm))

gplus <- sm / Nz
output <- paste0(args[3],".",args[4],".mds ",gplus)
write(output,file=paste0(args[1],"gplus_results.txt"),append=TRUE)
