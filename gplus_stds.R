#load accessory packages, set parameters
options(repos=structure(c(CRAN='http://cran.us.r-project.org')))
args <- commandArgs(trailingOnly=T)
cat("\n",args,as.character(Sys.time()),"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,doParallel,glmpca,irlba)

#location of data file
#args <- c('/users/ndyjack/Dist_Proj/tables/','/users/ndyjack/Dist_Proj/rdata/','malt_citeseq_2cl','glm_pca')
counts_loc <- paste0(args[1],args[3],"_expr.mtx") 
counts_mat <- readMM(counts_loc)
counts_mat <- t(as.matrix(counts_mat))
#counts_mat <- counts_mat[1:500,]
#counts_mat <- counts_mat[,-which(colSums(counts_mat)==0)]

#cell cluster identity vector
cluster_lab <- as.integer(readLines(paste0(args[1],args[3],"_labs.txt")))
#cluster_lab <- cluster_lab[1:500]
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

cat("\ndoing pca...\n")
if(args[4]=="std_pca"){
  pc_result <- as.matrix(prcomp(counts_mat)$x[,1:50])
}else if(args[4]=="trn_pca"){
  pc_result <- as.matrix(irlba(A=counts_mat,nv=50)$u)
}else if(args[4]=="glm_pca"){
  pc_result <- as.matrix(glmpca(Y=t(counts_mat),L=50,fam="mult",ctl = list(maxIter =100, eps = 1e-03))$factors)
}

dist_vec <- as.matrix(dist(pc_result))
dist_vec <- dist_vec[upper.tri(dist_vec)]

cat("\nestimating g+...\n")
n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
dist_inter <- sort(dist_vec[ind_inter])
dist_intra  <- sort(dist_vec[ind_intra])
rm(dist_vec) #memory cleanup
orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
gplus_samp <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nz_samp
#make this more general at some point, just set manually for now
#sm <- mclapply(1:length(ind_vec), function(i) {
#  tmp <- sum(sapply(ind_vec[[i]], function(j) sum(dist_vec[j] > dist_vec[ind_inter]))) 
  #isum(dist_vec[ind_intra[i]] > dist_vec[ind_inter])
#  return(tmp)
#}, mc.cores=30,mc.preschedule=F,mc.silent=F)
#sm <- sum(unlist(sm))
cat("\nfinished g+...\n")
#write output
#gplus <- sm/Nz
output <- paste0(args[3],".",args[4],".full ",gplus_samp)
cat("\n writing output to ",paste0(args[1],args[5]),"\n")
write(output,file=paste0(args[1],args[5]),append=TRUE)
