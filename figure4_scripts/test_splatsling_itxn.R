#mtold_conda_R/3.5.x
options("width"=180,warn=-1)
library(pacman)
p_load(DropletUtils,heatmap3,viridis,readxl,Matrix,scran,splatter,mbkmeans,scater,scry,purrr,scDesign2,SPARSim,irlba,slingshot)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colfun <- gg_color_hue(5)

cex_min <- 0.5
cex_max <- 2.5
#output = output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start)
map_cex <- function(x,y1,y2) {
  cex_min + ((cex_max - cex_min) / (y2 - y1) * (x - y1 ) )
}

plotme <- function(lab){

  labels <- dat@colData@listData$Group
  labcols <- sapply(as.numeric(labels), function(x) colfun[x])
  rna_counts <- as.matrix(dat@assays@data$counts)
  libsizes <- colSums(rna_counts)
  lim <- c(min(libsizes),max(libsizes))
  ptsizes <- sapply(libsizes, function(z) map_cex(z,lim[1],lim[2]))

  drop <- which(rowSums(rna_counts)==0)
  if(length(drop)>0){ rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),] }
  sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1),counts=rna_counts))
  fit <- modelGeneVar(sce)
  hvg_1k <- getTopHVGs(fit, n=1000)

  rna_counts <- t(rna_counts[hvg_1k,])
  tmp_dist <- dist(rna_counts,method = "manhattan")
  tmp_mds <- as.matrix(cmdscale(d=tmp_dist,k = 50)) 
  traj_lab <- as.matrix(as.numeric(labels))

  use <- which(!is.na(traj_lab[,1]))
  slingshot_res <- slingshot(data=tmp_mds[use,],clusterLabels=traj_lab[use,1],
    start.clus=min(traj_lab[use,1]),end.clus=max(traj_lab[use,1]),smoother='loess')
  lineages <- getLineages(slingshot_res)
  curves <- getCurves(lineages,reassign = TRUE,thresh=0.00001,reweight=TRUE,stretch=2,maxit=15)
  cors <- apply(slingPseudotime(curves),2,function(x) cor(x,traj_lab[use,1],use='complete.obs',method='spearman'))

  plot(x= tmp_mds[,1],y= tmp_mds[,2],main=lab, xlab="",ylab="",xaxt="n", yaxt="n",col=labcols,pch=16,cex=ptsizes)
  lines(curves)
  lines(lineages)
  mtext(side=1,text=formatC(cors,format='f',digits=4),line=0.5)
  mtext(side=2,paste(lineages@lineages$Lineage1,collapse=","),line=0.1)

}

#############################
#######simulate with splatter
#############################

#fixed
# k=5
# n=1000
# b=balanced


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
nC <- 5000
#k, number of celltypes
k <- 5
#cell type design
gPs <- rep(0.2,5)
#trajectory
t <- c(0,1,2,3,4)
#steps for continuous trajectory
s <- 1000

#easy medium and hard for lS, dP, and lF
lS_vec <- c(10.0,9.5,9.0)
dP_vec <- c(0.20,0.10,0.05)
#seq(0.05,0.45,by=0.05) #c(0.20,0.10,0.05)
lF_vec <- c(0.60,0.50,0.40)

#initial parameters
params <- newSplatParams()
params <- setParam(params, "nGenes" , nG)
params <- setParam(params, "seed", seed)
params <- setParam(params, "batchCells" , nC)
#params <- setParam(params, "steps" , s)

pdf("/users/ndyjack/Dist_Proj/images/splatsim_itxn_v8.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(2,2,2,1))

#vary lS, dP, and lF
params <- setParam(params, "de.facLoc", lF_vec[2])
for(i in 1:3){
  j <- l <- i
#  for(j in 1){
#    for(l in 1:length(dP_vec)){
      print(l)
      params <- setParam(params, "de.facLoc", lF_vec[i])
      params <- setParam(params, "lib.loc", lS_vec[j])
      dat <- splatSimulatePaths(params, group.prob = gPs, de.prob = dP_vec[l], path.from=t, verbose = FALSE,path.nSteps=s)
      tmplab <- paste0('lS',lS_vec[j],' lF',lF_vec[i], ' dP',dP_vec[l])
      plotme(lab=tmplab)
#    }
#  }
}

dev.off()
