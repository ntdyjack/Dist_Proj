#load accessory functions, set parameters
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,RcppCNPy,slingshot,reticulate)
np <-import("numpy")


#color label function for clusters
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colfun <- gg_color_hue(5)

#library size to point size function / labels
cex_min <- 0.5
cex_max <- 2.5
map_cex <- function(x,y1,y2) {
  cex_min + ((cex_max - cex_min) / (y2 - y1) * (x - y1 ) )
}
matrix_dir <- "/fastscratch/myscratch/ndyjack/simulations/"
mtx_files <- list.files(matrix_dir)
mtx_files <- grep('_t1_',mtx_files,value=T)
pt_cexs <- lapply(mtx_files, function(x) {
  dat <- readMM(paste0(matrix_dir,x))
  libsizes <- colSums(dat)
  lim <- c(min(libsizes),max(libsizes))
  ptsizes <- sapply(libsizes, function(z) map_cex(z,lim[1],lim[2]))
  return(ptsizes)
})

#3x3 plot, trajectory 1, easy medium hard, L1, LLR, HEL
labels_dir <- "/fastscratch/myscratch/ndyjack/labels/"
lab_files <- list.files(labels_dir)
lab_files <- grep("_t1_", lab_files,value=T)
pt_labs <- lapply(lab_files, function(x) {
  dat <- readLines(paste0(labels_dir,x))
  dat <- as.matrix(as.integer(as.factor(dat)))
  return(dat)
})
pt_cols <- lapply(pt_labs, function(x) sapply(x, function(y) colfun[y]))

dists_vec <- c('L1'=0, 'HEL'=6, 'LLR'=10)
pickles_dir <- "/fastscratch/myscratch/ndyjack/pickles/"

pdf("/users/ndyjack/Dist_Proj/images/troubleshoot_splat_trajcor_d2.pdf",width=12,height=12)
par(mfrow=c(3,3),mar=c(2,2,2,1))
for(i in 1:3){
  for(j in 1:3){

    tmp_nm <- paste0(pickles_dir,'splatsim_t1_d',j,'_hvg1k.mtx.',dists_vec[i],'.npy')
    tmp_dist <- np$load(tmp_nm)
    tmp_mds <- as.matrix(cmdscale(d=tmp_dist,k = 2))
    traj_lab <- pt_labs[[j]]
    use <- which(!is.na(traj_lab[,1]))
    slingshot_res <- slingshot(data=tmp_mds[use,],clusterLabels=traj_lab[use,1],
      start.clus=min(traj_lab[use,1]),end.clus=max(traj_lab[use,1]),smoother='loess')
    lineages <- getLineages(slingshot_res)
    curves <- getCurves(lineages)
    cors <- apply(slingPseudotime(curves),2,function(x) cor(x,traj_lab[use,1],use='complete.obs',method='spearman'))

    plot(x= tmp_mds[,1],y= tmp_mds[,2],main=paste(names(dists_vec[i]),'d',j),
      xlab="",ylab="",xaxt="n", yaxt="n",col=pt_cols[[j]],pch=16,cex=pt_cexs[[j]])
    lines(curves)
    mtext(side=1,text=formatC(cors,format='f',digits=4),line=0.5)
  }
}

dev.off()
