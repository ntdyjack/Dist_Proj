#load accessory packages, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,as.character(Sys.time()),"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,glmpca,irlba)

#location of data file
#args <- c("/users/ndyjack/Dist_Proj/tables/test_datasets/","/users/ndyjack/Dist_Proj/rdata/","cellbenchCS2_rnaseq_5cl_full","glm_pca")
counts_loc <- paste0(args[1],args[3],"_expr.mtx") 
counts_mat <- readMM(counts_loc)
counts_mat <- t(as.matrix(counts_mat))
#summary(colSums(counts_mat))
#summary(rowSums(counts_mat))
#dim(counts_mat)
#cluster_lab <- as.integer(readLines(paste0(args[1],args[3],"_labs.txt")))


if(args[4] %in% c("std_pca","trn_pca","glm_pca")){
  cat("running... ", args[4],"\n")
  start_time <- Sys.time()
  if(args[4]=="std_pca"){
    pc_result <- as.matrix(prcomp(counts_mat)$x[,1:50])
  }else if(args[4]=="trn_pca"){
    pc_result <- as.matrix(irlba(A=counts_mat,nv=50)$u)
  }else if(args[4]=="glm_pca"){
    pc_result <- as.matrix(glmpca(Y=t(counts_mat),L=50,fam="mult",ctl = list(maxIter =100, eps = 1e-03))$factors)
  }
  dist_tmp <- as.matrix(dist(pc_result))
  time_tmp <- difftime(time2=start_time,time1=Sys.time(),units="secs")
  dist_tmp <- dist_tmp[upper.tri(dist_tmp)]
  output_file <- paste0(args[2],args[3],".",args[4],".rda") 
  cat("writing results to... ",output_file,"\n\n")
  saveRDS(list(dist=dist_tmp,time=time_tmp),output_file)
} 

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQRL2'=2, 'JSM/PSM'=3,
  'JSD/PSD'=4, 'MKL'=5, 'POISSON'=6, 'HELLINGER'=7,
  'BHAT_METR'=8, 'BHAT_DIST'=9, 'TVD'=10,'LLR'=11,
  'REV_MKL'=12, 'REV_POIS'=13 ,'UWLLR'=14, 'ITA_SAI'=16,
  'REV_ITA_SAI'=17, 'COS_DIST'=18, 'PRB_COS_DIST'=19, 'COS_SIM'=20,
  'PRB_COS_SIM'=21,'PL2'=28,'PSL2'=29
)
i <- as.numeric(args[4])
j <-  which(dists_vec==i)
#metrics 16 and 17 need a pseudo count
if(i %in% c(5,6,12,13,16,17)){
  counts_mat <- counts_mat + 1.0
}
cat("running... ", dists_vec[j] ,names(dists_vec[j]),"\n")
start_time <- Sys.time()
dist_tmp <- Rfgc::dist_matrixdf(counts_mat, i)
time_tmp <- difftime(time2=start_time,time1=Sys.time(),units="secs")
dist_tmp <- dist_tmp[upper.tri(dist_tmp)] 

output_file <- paste0(args[2],args[3],".",args[4],".rda")
cat("writing results to... ",output_file,"\n\n") 
saveRDS(list(dist=dist_tmp,time=time_tmp),output_file)
