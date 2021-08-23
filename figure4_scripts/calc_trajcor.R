#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,RcppCNPy,slingshot,reticulate)
np <-import("numpy")

#args <- c("/users/ndyjack/Dist_Proj/pickles/","/users/ndyjack/Dist_Proj/tables/test_datasets/","rnamixss_rnaseq_traj_hvg500","0","trajcor_results_July2020.txt")
#cell cluster identity vector
#get rid of hvg string
traj_lab <- read.table(paste0(args[2],sub("_[h|f].*","",args[3]),"_labs.txt"),header=T)

#read in cell distance vector
#dist_vec <- readRDS(paste0(args[1],args[3],".",args[4],".rda"))$'dist'
#dist_mat <- matrix(0,nrow=nrow(traj_lab),ncol=nrow(traj_lab))
#dist_mat[upper.tri(dist_mat)] <- dist_vec
#dist_mat <- t(dist_mat)
#dist_mat[upper.tri(dist_mat)] <- dist_vec
dist_mat <- np$load(paste0(args[1],args[3],".",args[4],".npy"))

#make sure it's doing the right thing
#counts_mat <- as.matrix(readMM(paste0(args[2],args[3],"_expr.mtx")))
#tmp_dist <- Rfgc::dist_matrixdf(t(counts_mat), as.numeric(args[4]))
#all.equal(dist_mat,tmp_dist)

cat("\ncalculating MDS...\n")
mds_mat <- as.matrix(cmdscale(d=dist_mat,k = 50))

cat("\ncalculating slingshot...\n")
cors_by_traj <- sapply(1:ncol(traj_lab), function(i) {
  use <- which(!is.na(traj_lab[,i]))
  slingshot_res <- slingshot(data=mds_mat[use,],clusterLabels=traj_lab[use,i],start.clus=min(traj_lab[use,i]),end.clus=max(traj_lab[use,i]),smoother='loess')
  lineages <- getLineages(slingshot_res)
  curves <- getCurves(lineages,reassign = TRUE,thresh=0.00001,reweight=TRUE,stretch=2,maxit=15)
  pseudotime_cors <- apply(slingPseudotime(curves),2,function(x) cor(x,traj_lab[use,i],use='complete.obs',method='spearman'))
  return(max(abs(pseudotime_cors)))
})

output <- data.frame(pearson=cors_by_traj,trajectory=colnames(traj_lab))
outfile <- paste0(sub(x=args[2],pattern="test_datasets","trajcor_new"),args[3],".",args[4],".txt")
cat("\n writing output to ",outfile,"\n")
write.table(output,file=outfile,quote=F,row.names=F,col.names=T)
