#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,RcppCNPy,slingshot,reticulate)
np <-import("numpy")

#args <- c("/users/ndyjack/Dist_proj/pickles/","/users/ndyjack/Dist_Proj/tables/test_datasets/","rnamixss_rnaseq_traj_hvg500","0","trajcor_results_July2020.txt")
#args <- c("/fastscratch/myscratch/ndyjack/pickles/","/users/ndyjack/Dist_Proj/tables/test_datasets/","splatsim_t1_d2_hvg1k.mtx","0")
#cell cluster identity vector
#get rid of hvg string
#traj_lab <- readLines(paste0(args[2],sub("_[h|f].*","",args[3]),"_labs.txt"),header=T)
traj_lab <- as.matrix(as.integer(as.factor(readLines(paste0(sub('pickles','labels',args[1]),sub('hvg1k.mtx','labs.txt',args[3]))))))
if(!grepl('t1',args[3])){
  tmp <- matrix(nrow=nrow(traj_lab),ncol=2,NA)
  t1 <- which(traj_lab[,1] %in% c(1,2,3,4))
  t2 <- which(traj_lab[,1] %in% c(1,2,3,5))
  tmp[t1,1] <- traj_lab[t1,1]
  tmp[t2,2] <- traj_lab[t2,1]
  traj_lab <- tmp 
}

#read in cell distance vector
dist_mat <- np$load(paste0(args[1],args[3],".",args[4],".npy"))

cat("\ncalculating MDS...\n")
mds_mat <- as.matrix(cmdscale(d=dist_mat,k = 50))

cat("\ncalculating slingshot...\n")
cors_by_traj <- sapply(1:ncol(traj_lab), function(i) {
  use <- which(!is.na(traj_lab[,i]))
  slingshot_res <- slingshot(data=mds_mat[use,],clusterLabels=traj_lab[use,i],start.clus=min(traj_lab[use,i]),end.clus=max(traj_lab[use,i]),smoother='loess')
  lineages <- getLineages(slingshot_res)
  curves <- getCurves(lineages,reassign = TRUE,thresh=0.0001,reweight=TRUE,stretch=2,maxit=15)
  pseudotime_cors <- apply(slingPseudotime(curves),2,function(x) cor(x,traj_lab[use,i],use='complete.obs',method='spearman'))
  return(max(abs(pseudotime_cors)))
})

output <- data.frame(pearson=cors_by_traj,trajectory=1:ncol(traj_lab))
outfile <- paste0(sub(x=args[2],pattern="test_datasets","trajcor_new"),args[3],".",args[4],".txt")
cat("\n writing output to ",outfile,"\n")
write.table(output,file=outfile,quote=F,row.names=F,col.names=T)
