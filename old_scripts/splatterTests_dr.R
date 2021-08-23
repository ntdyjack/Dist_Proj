#load accessory functions
#args <- as.integer(commandArgs(trailingOnly=TRUE))
source("/users/ndyjack/Dist_Proj/scripts/functions.R")
library(splatter)
library(cluster)
library(PoiClaClu)

#set up parallel computing
library(doParallel)
cl <- makeCluster(3,type='FORK')
registerDoParallel(cl)

#splatter simulation paramters 
seed <- 1472
#total number of cells
nC <- 1000
#probability of being in each group
gP <- c(0.5,0.5)
#probability of difexp genes
dP <- 0.1
#number of genes
nG <- 20000
#dropout shape (lower values mean more dropout of dO is positive)
dS <- -100
#library size factor
lS <- 10

params <- newSplatParams()
params <- setParam(params, "nGenes" , nG)
params <- setParam(params, "seed", seed)
params <- setParam(params, "batchCells" , nC)
params <- setParam(params, "group.prob", gP)
params <- setParam(params, "de.prob", dP)
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.shape",dS)
params <- setParam(params, "lib.loc", lS)

#log fold change vector
lF_vec <- c(0,4,8) 
#c(args,args + 0.5) #seq(args,args+.9,length.out=2)

#dropout vectorr
dO_vec <- c(-700,-2,4) 

#parallelize over the dropout vector, since there's 10x as many of them
simulation_results <- foreach(i=1:length(dO_vec)) %dopar% {
#simulation_results <- lapply(lS_vec, function(lS) { 
  lapply(1:length(lF_vec), function(j) {
    dO <- dO_vec[i]
    lF <- lF_vec[j]
    #simulate data
    params <- setParam(params, "dropout.mid", dO)
    params <- setParam(params, "de.facLoc", lF)
    dat <-  splatSimulate(params,method="groups",verbose=F)
    group_clust <- dat@colData@listData$Group
    pct_dropout <- formatC( 100*(sum(dat@assays@data$Dropout) / (nC*nG)),
      format='f',digits=1)
    pct_zero <- formatC( 100*(sum(dat@assays@data$counts==0) / (nC*nG)),format='f',digits=1)
    dat <- dat@assays@data$counts
    info_lab <- paste0(nC, " cells, ",nG," genes, ~",median(colSums(dat))," UMIs, ",
      dP*100, "% DEGs, ",lF, "~LFC, ", pct_zero, "% Zeros, ",pct_dropout, "% Dropout" , dO)

    if(any(colSums(dat)==0)){
      group_clust <- group_clust[-which(colSums(dat)==0)]
      dat <- dat[,-which(colSums(dat)==0)]
    }

    cat(info_lab)
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

    return(list(euc=dist_euc,lrm=dist_lrm,poi=dist_poi,id=group_clust,info=info_lab))
  })
}
output_file <- paste0("/users/ndyjack/Dist_Proj/rdata/splatterTests.dr.rda")
saveRDS(simulation_results,output_file)
