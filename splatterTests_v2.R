#load accessory functions
args <- as.integer(commandArgs(trailingOnly=TRUE))
source("/users/ndyjack/Dist_Proj/scripts/functions.R")
library(splatter)
library(cluster)
library(PoiClaClu)

#set up parallel computing
library(doParallel)
cl <- makeCluster(10,type='FORK')
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
lF_vec <- c(args,args + 0.5) #seq(args,args+.9,length.out=2)

#dropout vectorr
dO_vec <- c(-700,-300,-125,-75,-50,-35,-25,-15,-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4) 
  #c(seq(-700,-200,length.out=5),seq(-180,-70,length.out=5),seq(-65,-30,length.out=10),
  #seq(-28,-10,length.out=15),seq(-9,-5,length.out=10),seq(-4,1,length.out=50),seq(1,4,length.out=5))
#dO_vec <- log(c(seq(1e-300,1e-200,length.out=24),seq(1e-199,1e-100,length.out=23),
#    seq(1e-99,1e-50,length.out=22),seq(1e-49,1e-1,length.out=21),seq(1e0,60,length.out=10))) 
#seq(-500,5,length.out=100) 
#log(c(seq(),seq())) 
#(1/2) ^ seq(-10,10,length.out=100)

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

#   cat(info_lab)
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
    return(list(sil=sil,info=info_lab,times=c(time_euc,time_lrm,time_poi)))
  })
}
output_file <- paste0("/users/ndyjack/Dist_Proj/rdata/splatterTests.lF.",args,".rda")
saveRDS(simulation_results,output_file)
