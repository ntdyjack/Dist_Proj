#options("width"=160)
args <- as.integer(commandArgs(trailingOnly=TRUE))
source("/users/ndyjack/Dist_Proj/scripts/functions.R")
library(PoiClaClu)
library(DropletUtils)

#load in null data
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19") 
counts <- as.matrix(counts(sce))
#drop rows with zero expression
counts <- counts[-which(rowSums(counts)==0),]
#drop cells with CLEC9A expression (these are CD14+ DC's, mess up homogenity)
counts <- counts[,-which(counts["ENSG00000197992",]>0)]
#compute variances for each gene
gene_vars <- apply(counts,1,var)

#set up parallel computing
library(doParallel)
cl <- makeCluster(5,type='FORK')
registerDoParallel(cl)


#gene breakdown for parallelization
n_genes <- nrow(counts)
n_jobs <- 15
bin_size <- n_genes/n_jobs
n_pars <- 20
nG_vec <- round(seq(bin_size*args + 1, bin_size*(1+args),length=n_pars))
#lapply(0:(n_jobs-1), function(args) round(seq(bin_size*args + 1, bin_size*(1+args),length=20)))

simulation_results <- foreach(i=1:length(nG_vec)) %dopar% {
#simulation_results <- lapply(lS_vec, function(lS) {
    dat <- counts[order(gene_vars,decreasing=T)[1:i],]
 
    #test distances
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
      euc = silhouette(x=as.integer(group_clust),dist=dist_euc)[,3],
      lrm = silhouette(x=as.integer(group_clust),dist=dist_lrm)[,3],
      poi = silhouette(x=as.integer(group_clust),dist=dist_poi)[,3]
    )
#})
    return(list(sil=sil,info=info_lab,times=c(time_euc,time_lrm,time_poi)),nGenes=i)
}
output_file <- paste0("/users/ndyjack/Dist_Proj/rdata/cd14.nG.",args,".rda")
saveRDS(simulation_results,output_file)
