#mtold_conda_R/3.5.x
options("width"=180,warn=-1)
library(pacman)
p_load(DropletUtils,heatmap3,viridis,readxl,Matrix,scran,splatter,mbkmeans,scater,scry,purrr,scDesign2,SPARSim,irlba)
output_dir <- "/fastscratch/myscratch/ndyjack/simulations/splatsim"
#bname <- '/fastscratch/myscratch/ndyjack/simulations/splatsim'

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cex_min <- 0.5
cex_max <- 2.5
#output = output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start)
map_cex <- function(x,y1,y2) {
  cex_min + ((cex_max - cex_min) / (y2 - y1) * (x - y1 ) )
}

plotme <- function(x=NULL,y=NULL){
  labels <- dat@colData@listData$Group
  labcols <- sapply(as.numeric(labels), function(x) colfun[x])
  outlier <- rowData(dat)[,'OutlierFactor']!=1
  rna_counts <- as.matrix(dat@assays@data$counts)
  #lim <- c(min(rna_counts),max(rna_counts))
  libsizes <- colSums(rna_counts)
  lim <- c(min(libsizes),max(libsizes))
  ptsizes <- sapply(libsizes, function(z) map_cex(z,lim[1],lim[2]))
  #rm(dat) #cleanup
  drop <- which(rowSums(rna_counts)==0)
  if(length(drop)>0){ rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),] }
  sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1),counts=rna_counts))
  #rm(rna_counts) #cleanup
  fit <- modelGeneVar(sce)
  hvg_1k <- getTopHVGs(fit, n=1000)
  sce <- computeSumFactors(sce, min.mean=0.1)
  sce <- logNormCounts(sce)
  sce <- scater::runPCA(sce, ncomponents = 50,subset_row=hvg_1k,scale = TRUE, name='hvg')
  #sce <- scater::runPCA(sce, ncomponents = 50,ntop = nrow(rna_counts),scale = TRUE, name='ful')
  plot(x= reducedDims(sce)[[1]][,1],y= reducedDims(sce)[[1]][,2],main=x,
    xlab="",ylab="",xaxt="n", yaxt="n",
    col=labcols,pch=ifelse(outlier,13,16),cex=ptsizes)
  mtext(text=y,side=2,line=0.0)
  #plot(x= reducedDims(sce)[[2]][,1],y= reducedDims(sce)[[2]][,2],main="full",xlab="",ylab="",xaxt="n", yaxt="n",col=labcols,pch=16)
  rna_counts <- Matrix(counts(sce)[hvg_1k,],sparse=T)
  #writeMM(obj = rna_counts, file=paste0(oname,"_full.mtx"))
  writeMM(obj = rna_counts, file=paste0(oname,"_hvg1k.mtx"))
  writeLines(as.character(labels),paste0(sub('simulations','labels',oname),"_labs.txt"))
  #rm(sce) 
}

#############################
#######simulate with splatter
#############################

#explore
# k (number clusters, 2,5,10)
# n (number cells, 1k, 5k, 10k)
# b, % celltypes (balanced, medium, unbalanced)


#roll into easy medium and hard
# l, library size (default, cut in half, double ~ make them around 100 to low 1k)
# p, % of differentially expressed genes (5, 10, 20 %?)
# f, fold change (default, x2, and x4 of that)


#############################
#######splatter-simulated data
#############################

#splatter simulation paramters
seed <- 1472
#number of genes
nG <- 20000

#n, total number of cells
nC_vec <- c(1000,5000,10000)

#k, number of celltypes
k_vec <- c(2,5,10)

#b, balance of celltypes
#stable design for k=2,5,10
gPs_list <- list(
  rep(0.5,2),
  rep(0.2,5),
  rep(0.1,10)
)
#variable design for k=5
gPv_list <- list(
  c(0.20,0.20,0.20,0.20,0.20), #most balanced
  c(0.15,0.35,0.25,0.15,0.10), #medium balanced
  c(0.05,0.55,0.25,0.05,0.10) #less balanced
)

#t, trajectories?
t_list <- list(
  c(0,1,2,3,4),
  c(0,1,2,3,3),
  c(0,1,2,3,3)
)
gPt_list <- list(
  c(0.20,0.20,0.20,0.20,0.20),
  c(0.20,0.20,0.20,0.20,0.20),
  c(0.20,0.20,0.20,0.35,0.05)
)

#easy medium and hard for lS, dP, and lF
lS_vec <- c(10.0,9.5,9.0)
dP_vec <- c(0.20,0.10,0.05)
lF_vec <- c(0.30,0.20,0.10)
lFt_vec <- c(1.0,0.5,0.25) #trajectory settings need to be easier
set_vec <- c('Easy','Medium','Hard')

#steps for continuous trajectory
s <- 2000

#initial parameters
params <- newSplatParams()
params <- setParam(params, "nGenes" , nG)
params <- setParam(params, "seed", seed)

#bname <- paste0(output_dir,'splatsim_vary')
#tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[i],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[2])

pdf("/users/ndyjack/Dist_Proj/images/splattersim_pcaplots_v7.pdf",width=12,height=12)
par(mfrow=c(3,3),mar=c(2,2,2,1))

#vary n with fixed k=5 and b=most
#coloring function
colfun <- gg_color_hue(k_vec[2])
for(i in 1:3){
#  print(i)
  for(j in 1:3){
    print(paste0(i,',',j))
    params <- setParam(params, "batchCells" , nC_vec[i])
    params <- setParam(params, "lib.loc", lS_vec[j])
    params <- setParam(params, "de.facLoc", lF_vec[j])
    dat <- splatSimulateGroups(params, group.prob = gPs_list[[2]], de.prob = rep(dP_vec[j],k_vec[2]),verbose = FALSE)
    oname <- paste0(output_dir,'_n',i,'_d',j)
    #paste0(bname[2],nC_vec[2],bname[3],k_vec[i],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[2])
    plotme(x=set_vec[j],y=paste0('n',nC_vec[i]))
  }
}

#vary k with fixed n=5000 and b=most
params <- setParam(params, "batchCells" , nC_vec[2])
for(i in 1:3){
  #print(i)
  for(j in 1:3){
    print(paste0(i,',',j))
    colfun <- gg_color_hue(k_vec[i])
    params <- setParam(params, "lib.loc", lS_vec[j])
    params <- setParam(params, "de.facLoc", lF_vec[j])
    dat <- splatSimulateGroups(params, group.prob = gPs_list[[i]], de.prob = rep(dP_vec[j],k_vec[i]),verbose = FALSE)
    oname <- paste0(output_dir,'_k',i,'_d',j)
    plotme(x=set_vec[j],y=paste0('k',k_vec[i]))
  }
}

#vary b with fixed n=5000 and k=5
colfun <- gg_color_hue(k_vec[2])
params <- setParam(params, "batchCells" , nC_vec[2])
for(i in 1:3){
  for(j in 1:3){
    print(paste0(i,',',j))
    params <- setParam(params, "lib.loc", lS_vec[j])
    params <- setParam(params, "de.facLoc", lF_vec[j])
    dat <- splatSimulateGroups(params, group.prob = gPv_list[[i]], de.prob = rep(dP_vec[j],k_vec[2]),verbose = FALSE)
    oname <- paste0(output_dir,'_b',i,'_d',j)
    plotme(x=set_vec[j],y=paste0('b',i))
  }
}

#vary t(trajectory) with fixed n=5000
params <- setParam(params, "batchCells" , nC_vec[2])
for(i in 1:3){
  for(j in 1:3){
    print(paste0(i,',',j))
    params <- setParam(params, "lib.loc", lS_vec[j])
    params <- setParam(params, "de.facLoc", lFt_vec[j])
    dat <- splatSimulatePaths(params, group.prob = gPt_list[[i]], path.from=t_list[[i]], de.prob = dP_vec[j], verbose = FALSE, path.nSteps=s)
    oname <- paste0(output_dir,'_t',i,'_d',j)
    plotme(x=set_vec[j],y=paste0('t',i))
  }
}

dev.off()

#for(i in 1:3){
#  colfun <- gg_color_hue(k_vec[i])
#  tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[i],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[2])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  dat <- splatSimulateGroups(params, group.prob = gPs_list[[i]], de.prob = rep(dP_vec[2],k_vec[i]),verbose = FALSE)
# plotme()
#}

#dev.off()

#colfun <- gg_color_hue(k_vec[2])
#vary n, number cells
#for(i in 1:3){
#  tname <- paste0(bname[2],nC_vec[i],bname[3],k_vec[2],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[2])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  params <- setParam(params, "batchCells" , nC_vec[i])
#  dat <- splatSimulateGroups(params, group.prob = gPs_list[[2]], de.prob = rep(dP_vec[2],k_vec[2]),verbose = FALSE)
#  plotme()
#}

#params <- setParam(params, "batchCells" , nC_vec[2])
#vary b, balance of the design
#for(i in 1:3){
#  tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[2],bname[4],i,bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[2])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  dat <- splatSimulateGroups(params, group.prob = gPv_list[[i]], de.prob = rep(dP_vec[2],k_vec[2]),verbose = FALSE)
#  plotme()
#}

#vary l, library size
#for(i in 1:3){
#  tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[2],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[i],bname[7],lF_vec[2])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  params <- setParam(params, "lib.loc" , lS_vec[i])
#  dat <- splatSimulateGroups(params, group.prob = gPs_list[[2]], de.prob = rep(dP_vec[2],k_vec[2]),verbose = FALSE)
#  plotme()
#}

#params <- setParam(params, "lib.loc" , lS_vec[2])
#vary p, percentage of differentially expressed genes
#for(i in 1:3){
#  tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[2],bname[4],'1',bname[5],dP_vec[i],bname[6],lS_vec[2],bname[7],lF_vec[2])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  dat <- splatSimulateGroups(params, group.prob = gPs_list[[2]], de.prob = rep(dP_vec[i],k_vec[2]),verbose = FALSE)
#  plotme()
#}

#vary f, fold-change
#for(i in 1:3){
#  tname <- paste0(bname[2],nC_vec[2],bname[3],k_vec[2],bname[4],'1',bname[5],dP_vec[2],bname[6],lS_vec[2],bname[7],lF_vec[i])
#  oname <- paste0(output_dir,bname[1],tname)
#  tname <- gsub("_"," ",tname)
#  params <- setParam(params, "de.facLoc" , lF_vec[i])
#  dat <- splatSimulateGroups(params, group.prob = gPs_list[[2]], de.prob = rep(dP_vec[2],k_vec[2]),verbose = FALSE)
#  plotme()
#}

#dev.off()

