##
#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,heatmap3,beeswarm,doParallel)

#parallel backend
options(cores=30)
registerDoParallel(cores=30)

dists_vec <- c(
  'L1','L2','SQRL2','JSM/PSM','JSD/PSD','MKL',
  'POISSON','HELLINGER','BHAT_METR','BHAT_DIST',
  'TVD','LLR','EMD','REV_MKL','REV_POISSON','UWLLR','OLLR'
)
idx_r <- 1:length(dists_vec)
idx_r <- idx_r[-which(dists_vec=="EMD")]
idx_c <- idx_r - 1

#location of data files
data_dir <- "/users/ndyjack/Dist_Proj/rdata/"
data_names <- list.files(data_dir)
data_names <- unique(sapply(data_names, function(x) strsplit(x,'\\.')[[1]][1]))
rdata_names <- lapply(data_names, function(x) paste(x,idx_c,"rda",sep="."))
#location of identity files
ident_dir <- "/users/ndyjack/Dist_Proj/tables/"
ident_vec <- sapply(data_names, function(x) paste0(x,"_labs.txt")) 
#location of output files
output_dir <- "/users/ndyjack/Dist_Proj/images/"





#loop over each dataset
output_metrics <- sapply(1:length(ident_vec), function(i) {
  #read in identity vector
  cluster_lab <- as.numeric(readLines(paste0(ident_dir,ident_vec[i])))
  Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
  Nz <- (Nt*(Nt-1))/2
  #coerce to paired vector (IE, upper triangular cluster idenity matrix)
  ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
  ind_vec <- apply(ind_vec,1,function(x) {
    tmp <- rep(0,length(x))
    tmp[which(x)] <- cluster_lab[which(x)]
    return(tmp)
  })
  ind_vec <- ind_vec[upper.tri(ind_vec)]
  ind_intra <- which(ind_vec>0)
  ind_inter <- which(ind_vec==0) 
  #loop over each metric
  output_tmp <- sapply(1:5, function(j) {
    cat(i,j)
    #read in distance vector
    dist_vec <- readRDS(paste0(data_dir,rdata_names[[i]][j]))[[1]]
    #s- number of times intra-clsuter distances are strictly greater than inter-cluster-distances
    sm <- foreach(k=ind_intra,.combine=c) %dopar% { sum(dist_vec[k] > dist_vec[ind_inter]) }
    #G-plus index
    return(sum(sm) / Nz)
  })
  return(output_tmp)
})
