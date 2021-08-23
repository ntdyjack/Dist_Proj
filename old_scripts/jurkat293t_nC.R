#options("width"=160)
args <- as.integer(commandArgs(trailingOnly=TRUE))
source("/users/ndyjack/Dist_Proj/scripts/functions.R")
library(PoiClaClu)
library(DropletUtils)
library(readxl)
library(cluster)
library(mclust)

#load in data
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/jurkat_293t_50_50_filtered_gene_bc_matrices/filtered_matrices_mex/hg19") 
counts <- as.matrix(counts(sce))
#drop rows with zero expression
counts <- counts[-which(rowSums(counts)==0),]
#make identity vector for jurkat vs 239t
group_clust <- read_excel("/users/ndyjack/Dist_Proj/tables/41467_2017_BFncomms14049_MOESM830_ESM.xlsx")
group_clust <- as.integer(as.factor(group_clust$`Species Call by Expression`))

#set up parallel computing
library(doParallel)
cl <- makeCluster(5,type='FORK')
registerDoParallel(cl)

n_cells <- ncol(counts)
n_jobs  <- 10
bin_size <- n_cells/n_jobs
n_pars <- 15
nC_mat <- matrix(round(seq(10,n_cells,length.out=n_pars * n_jobs)), ncol=n_pars, nrow=n_jobs)
if(nC_mat[n_jobs,n_pars] != n_cells){
    nC_mat <- nC_mat + (n_cells - nC_mat[n_jobs,n_pars])
}
nC_vec <- nC_mat[args+1,]
seed <- 1752

#lapply(0:(n_jobs-1), function(args) round(seq(bin_size*args + 1, bin_size*(1+args),length=20)))
#sapply(0:(n_jobs-1), function(args) round(seq(bin_size*args + 2, bin_size*(1+args),length=15)))
#compute variances for each gene
#gene_vars <- apply(counts,1,var)
#gene breakdown for parallelization
#n_genes <- nrow(counts)
#n_jobs <- 16
#bin_size <- n_genes/n_jobs
#n_pars <- 20
#nG_vec <- round(seq(bin_size*args + 1, bin_size*(1+args),length=n_pars))
#if(any(nG_vec==1)){ nG_vec[nG_vec==1] <- 2 }
#lapply(0:(n_jobs-1), function(args) round(seq(bin_size*args + 1, bin_size*(1+args),length=20)))

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
