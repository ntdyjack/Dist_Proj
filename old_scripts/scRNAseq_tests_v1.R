#options("width"=160)
args <- as.integer(commandArgs(trailingOnly=TRUE))
filedat <- read.csv("")
source("/users/ndyjack/Dist_Proj/scripts/functions.R")
library(PoiClaClu)
library(scRNAseq)
library(mclust)
library(cluster)

fun <- get(filedat$dat, mode = "function", envir = parent.frame(1)) 
sce <- fun()
counts <- as.matrix(assays(sce)[[1]])
counts <- counts[-which(rowSums(counts)==0),]
group_clust <- as.integer(as.factor(sce@colData[[filedat$outcome]]))

#set up parallel computing
library(doParallel)
cl <- makeCluster(5,type='FORK')
registerDoParallel(cl)

seed <- 1547

n_cells <- ncol(counts)
nC_vec <- round(seq(0.05,1.0,by=0.05) * n_cells) 
#round(seq(10,n_cells,length.out=n_pars))

#simulation_results <- foreach(i=1:5) %dopar% {
simulation_results <- foreach(i=1:length(nC_vec)) %dopar% {
#simulation_results <- lapply(lS_vec, function(lS) {
    #subset to relevant genes
    set.seed(seed)
    cells_use <- sample(1:n_cells,size=nC_vec[i])
    group_labs <- group_clust[cells_use] 
    dat <- counts[,cells_use]

    #multinomial likelihood ratio
    start_time <- Sys.time()
    dist_lrm <- compute_distance(dat,"lrm")
    time_lrm <- Sys.time() - start_time
    #euclidean distance (unscaled)
    start_time <- Sys.time()
    dist_euc <- compute_distance(dat,"euc")
    time_euc <- Sys.time() - start_time
    #Whitten Poisson Distance
    start_time <- Sys.time()
    dist_poi <- PoissonDistance(t(dat))$dd
    time_poi <- Sys.time() - start_time

    #silhouette scores
    sil <- data.frame(
      euc = silhouette(x=group_labs,dist=dist_euc)[,3],
      lrm = silhouette(x=group_labs,dist=dist_lrm)[,3],
      poi = silhouette(x=group_labs,dist=dist_poi)[,3]
    )

    #kmedoids cluster
    clust <- data.frame(
      euc = pam(dist_euc,k=2)$cluster, 
      lrm = pam(dist_lrm,k=2)$cluster,
      poi = pam(dist_poi,k=2)$cluster
    )

    #adjusted rand indexfrom the kmedoids cluster
    ari <- apply(clust,2,function(x) adjustedRandIndex(x,group_labs))
#})
    return(list(sil=sil,clust=clust,ari=ari,times=c(time_euc,time_lrm,time_poi),nCells=nC_vec[i],labs=group_labs))
}
output_file <- paste0("/users/ndyjack/Dist_Proj/rdata/jurkat293t.nC.",args,".rda")
saveRDS(simulation_results,output_file)
