args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet,RcppCNPy,reticulate,Matrix,cluster,plotrix,TeachingDemos,wordcloud)
np <- import("numpy")

#Figure 1A, MDS plots for LLR and L2, (real data hvg1k), colored by library size and PAM(k=2) (2x2)
#Figure 1B, MDS plots for LLR and L2, (sim data, hvg1k), colored by library sized and PAM(k-2) (2x2)
#Figure 1C, G(k) for all distances (real data hvg1k)
#Figure 1D, G(k) for all distances (sim data hvg1k)
#Figure 1 layout
#legends (libsize, cluster, distances)
#A, C
#B, D

set.seed(15123)
k <- 2

#load relevant distance matrices
pickle_dir  <-  "/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_hvg1k"
rna_lr_dist <- as.dist(np$load(paste0(pickle_dir,'.10.npy')))
rna_l2_dist <- as.dist(np$load(paste0(pickle_dir,'.1.npy')))

pickle_dir  <-  "/fastscratch/myscratch/ndyjack/pickles/splatsim_1cl_hvg1k"
sim_lr_dist <- as.dist(np$load(paste0(pickle_dir,'.10.npy')))
sim_l2_dist <- as.dist(np$load(paste0(pickle_dir,'.1.npy')))

#induce PAM for k=2, get labels
cl_rl2 <- pam(x=rna_l2_dist,k=k,diss=T)$clustering
cl_rlr <- pam(x=rna_lr_dist,k=k,diss=T)$clustering
cl_sl2 <- pam(x=sim_l2_dist,k=k,diss=T)$clustering
cl_slr <- pam(x=sim_lr_dist,k=k,diss=T)$clustering

cl_col_rl2 <- ifelse(cl_rl2==1,'#ff000064','#0000ff64')
cl_col_rlr <- ifelse(cl_rlr==1,'#ff000064','#0000ff64')
cl_col_sl2 <- ifelse(cl_sl2==1,'#ff000064','#0000ff64')
cl_col_slr <- ifelse(cl_slr==1,'#ff000064','#0000ff64')

#calculate 2d MDS for each
mds_rl2 <- as.matrix(cmdscale(d=rna_l2_dist,k = 2))
mds_rlr <- as.matrix(cmdscale(d=rna_lr_dist,k = 2))
mds_sl2 <- as.matrix(cmdscale(d=sim_l2_dist,k = 2))
mds_slr <- as.matrix(cmdscale(d=sim_lr_dist,k = 2))

#mds plot limits
xlim_rl2 <- c(min(mds_rl2[,1]),max(mds_rl2[,1]))
ylim_rl2 <- c(min(mds_rl2[,2]),max(mds_rl2[,2]))

xlim_rlr <- c(min(mds_rlr[,1]),max(mds_rlr[,1]))
ylim_rlr <- c(min(mds_rlr[,2]),max(mds_rlr[,2]))

xlim_sl2 <- c(min(mds_sl2[,1]),max(mds_sl2[,1]))
ylim_sl2 <- c(min(mds_sl2[,2]),max(mds_sl2[,2]))

xlim_slr <- c(min(mds_slr[,1]),max(mds_slr[,1]))
ylim_slr <- c(min(mds_slr[,2]),max(mds_slr[,2]))

#load original expression matrices to calculate library sizes factors
count_dir  <-  "/users/ndyjack/Dist_Proj/tables/test_datasets/only293t_rnaseq_1cl_full_expr.mtx"
count_mat <- readMM(count_dir)
cs_rna <- colSums(count_mat)
zlim_rna <- c(min(cs_rna),max(cs_rna))

count_dir  <-  "/users/ndyjack/Dist_Proj/tables/test_datasets/splatsim_1cl_full_expr.mtx"
count_mat <- readMM(count_dir)
cs_sim <- colSums(count_mat)
zlim_sim <- c(min(cs_sim),max(cs_sim))

#color palettes
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- heat.colors(nlevs-1) #viridis(nlevs-1)

col_pal_rna <- colorNumeric(col_pal, domain=zlim_rna, na.color = "#808080",alpha = FALSE,reverse = FALSE)
cs_col_rna <- col_pal_rna(cs_rna)
cs_col_rna <- paste0(cs_col_rna,'64')

col_pal_sim <- colorNumeric(col_pal, domain=zlim_sim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
cs_col_sim <- col_pal_rna(cs_sim)
cs_col_sim <- paste0(cs_col_sim,'64')

zleglocs <- c(0.01,0.99)
zleglabs <- c("min","max")

#gap statistics
#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)
l <- length(dists_vec)

palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))
pch_labs <- rep(c(15,16,17),5)

xlim_gs <- c(0.5,5.5)
xticks_gs <- 1:5

map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/hvg1k.'
gs_rna <- sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  return(tmp)  
})
gs_rna_scaled <- apply(gs_rna,2, function(x) map_fxn(x,0,1,min(x),max(x)))

results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/splat.hvg1k.'
gs_sim <- sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  return(tmp)
})
gs_sim_scaled <- apply(gs_sim, 2, function(x) map_fxn(x,0,1,min(x),max(x)))


pdf("/users/ndyjack/Dist_Proj/figures_final/figure_1.2A.pdf",width=4,height=4)
  plot.new()
  #color bar legend for library size
  par(new = "TRUE",plt = c(0.30,0.70,0.92,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  text(x=0.5,y=0.5,labels='Library Size', cex=0.8)
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.4,las=1,mgp=c(1.0, .4, 0))

  #labels for far left
  par(new = "TRUE",plt = c(0.00,0.06,0.01,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=rep(0.5,2),y=c(0.60,0.20),labels=c('scRNAseq','Splatter'),cex=1.2,srt=90)

  #Figure 1A.1 (topleft), real L2 colored by library size
  par(new = "TRUE",plt = c(0.15,0.50,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_rl2, ylim = ylim_rl2, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_rl2[,1],y=mds_rl2[,2],col=cs_col_rna,cex=0.8,pch=16)
  mtext(side=3,text="L2",cex=1.0,line=0.0,las=1)
  mtext(side=2,text="MDS2",cex=0.7,line=0.0,las=3)

  #Figure 1A.2 (topright), real LR colored by library size
  par(new = "TRUE",plt = c(0.55,0.90,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_rlr, ylim = ylim_rlr, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_rlr[,1],y=mds_rlr[,2],col=cs_col_rna,cex=0.8,pch=16)
  mtext(side=3,text="LLR",cex=1.0,line=0.0,las=1)

  #Figure 1A.3 (bottomleft), sim L2 colored by library size
  par(new = "TRUE",plt = c(0.15,0.50,0.05,0.40),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_sl2, ylim = ylim_sl2, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_sl2[,1],y=mds_sl2[,2],col=cs_col_sim,cex=0.8,pch=16)
  mtext(side=2,text="MDS2",cex=0.7,line=0.0,las=3)
  mtext(side=1,text="MDS1",cex=0.7,line=0.0,las=1)

  #Figure 1A.4 (bottomright), sim LR colored by library size
  par(new = "TRUE",plt = c(0.55,0.90,0.05,0.40),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_slr, ylim = ylim_slr, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_slr[,1],y=mds_slr[,2],col=cs_col_sim,cex=0.8,pch=16)
  mtext(side=1,text="MDS1",cex=0.7,line=0.0,las=1,mgp=c(0, 0, 0))

dev.off()

pdf("/users/ndyjack/Dist_Proj/figures_final/figure_1.2B.pdf",width=4,height=4)
  plot.new()
  #plot for PAM clustering
  par(new = "TRUE",plt = c(0.30,0.70,0.92,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  points(x=c(0.2,0.7),y=c(0.70,0.70),pch=16,cex=1.5,col=c('blue','red'))
  text(x=c(0.2,0.7),y=c(0.20,0.20),labels=c('Cl1 (PAM)', 'Cl2 (PAM)'), cex=0.8)

  #Figure 1A.1 (topleft), real L2 colored by PAM
  par(new = "TRUE",plt = c(0.15,0.50,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_rl2, ylim = ylim_rl2, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_rl2[,1],y=mds_rl2[,2],col=cl_col_rl2,cex=0.8,pch=16)
  mtext(side=3,text="L2",cex=1.0,line=0.0,las=1)

  #Figure 1A.2 (topright), real LR colored by PAM
  par(new = "TRUE",plt = c(0.55,0.90,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_rlr, ylim = ylim_rlr, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_rlr[,1],y=mds_rlr[,2],col=cl_col_rlr,cex=0.8,pch=16)
  mtext(side=3,text="LLR",cex=1.0,line=0.0,las=1)

  #Figure 1A.3 (bottomleft), sim L2 colored by PAM
  par(new = "TRUE",plt = c(0.15,0.50,0.05,0.40),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_sl2, ylim = ylim_sl2, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_sl2[,1],y=mds_sl2[,2],col=cl_col_sl2,cex=0.8,pch=16)
  mtext(side=1,text="MDS1",cex=0.7,line=0.0,las=1)

  #Figure 1A.4 (bottomright), sim LR colored by PAM
  par(new = "TRUE",plt = c(0.55,0.90,0.05,0.40),las = 1,cex.axis = 1)
  plot.window(xlim = xlim_slr, ylim = ylim_slr, xaxs = "i",yaxs = "i")
  box()
  points(x=mds_slr[,1],y=mds_slr[,2],col=cl_col_slr,cex=0.8,pch=16)
  mtext(side=1,text="MDS1",cex=0.7,line=0.0,las=1,mgp=c(0, 0, 0))

dev.off()

  #sprd_xnmz1 <- rep(2,length(dists_vec))
  #sprd_xnmz2 <- sprd_xnmz1
  #sprd_xnmz2[seq(1,l,2)] <- 2.5
  #sprd_xnmz2[seq(2,l-1,2)] <- 1.5

#0000007b

set.seed(123)
pdf("/users/ndyjack/Dist_Proj/figures_final/figure_1.2CD.pdf",width=4,height=4)
  plot.new()
  #distance type legend
  par(new = "TRUE",plt = c(0.10,0.95,0.80,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.6,y.intersp=1.1,text.width=0.15,cex=0.7,
    pch=pch_labs,pt.cex=1.3,bty='n')

  #Figure C
  par(new = "TRUE",plt = c(0.15,0.95,0.45,0.80))
  plot.new()
  plot.window(xlim = xlim_gs, ylim = c(-0.05,1.05), xaxs = "i",yaxs = "i")
  for(i in 1:l){
    points(x=1:5,y=gs_rna_scaled[,i],type='b',pch=pch_labs[i],col=col_labs[i],cex=0.8)
  }
  axis(side=2,labels=zleglabs,at=zleglocs,cex.axis=0.7,las=2,mgp=c(3, .5, 0))
  mtext(side=2,text="Gap Statistic (scaled)",cex=0.7,line=2,las=3)
  labu <- which(!names(dists_vec) %in% c('LLR','ISD'))
  xlbu <- sample(c(2,3),length(labu),replace=T)
  tmp_coords <- wordlayout(x=xlbu, y=gs_rna_scaled[2,labu],
    word = names(dists_vec)[labu] ,xlim=c(1.00,5.00), ylim= c(-0.01,1.01), tstep=0.25, rstep=0.25)
  text(x=tmp_coords[,1],y=tmp_coords[,2],labels=names(dists_vec)[labu],cex=0.5)
  for(i in 1:length(labu)){
    lines(x=c(xlbu[i],tmp_coords[i,1]),y=c(gs_rna_scaled[xlbu[i],labu[i]],tmp_coords[i,2]),lwd=0.5,type='l',col='#0000007b')
  }

  #Figure D
  par(new = "TRUE",plt = c(0.15,0.95,0.05,0.40))
  plot.new()
  plot.window(xlim = xlim_gs, ylim = c(-0.05,1.05), xaxs = "i",yaxs = "i")
  for(i in 1:length(dists_vec)){
    points(x=1:5,y=gs_sim_scaled[,i],type='b',pch=pch_labs[i],col=col_labs[i],cex=0.8)
  }
  axis(side=1,labels=xticks_gs,at=xticks_gs,cex.axis=0.6,las=1,mgp=c(3, .15, 0))
  axis(side=2,labels=zleglabs,at=zleglocs,cex.axis=0.7,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text="Gap Statistic (scaled)",cex=0.7,line=2,las=3)
  labu <- which(!names(dists_vec) %in% c('LLR','TVD','ISD','RISD'))
  xlbu <- sample(c(2,3),length(labu),replace=T)
  tmp_coords <- wordlayout(x=xlbu, y=gs_sim_scaled[2,labu],
    word = names(dists_vec)[labu] ,xlim=c(1.00,5.00), ylim= c(-0.01,1.01), tstep=0.25, rstep=0.25)
  text(x=tmp_coords[,1],y=tmp_coords[,2],labels=names(dists_vec)[labu],cex=0.5)
  for(i in 1:length(labu)){
    lines(x=c(xlbu[i],tmp_coords[i,1]),y=c(gs_sim_scaled[xlbu[i],labu[i]],tmp_coords[i,2]),lwd=0.5,type='l',col='#0000007b')
  }
dev.off()


nsim <- 500
xlim <- ylim <- c(-3,3)
cl1 <- '#00ff0064'
cl2 <- '#0000ff64'
cl3 <- '#ff000064'
set.seed(123)
d1 <- data.frame(x=rnorm(n=nsim,sd=0.5),y=rnorm(n=nsim,sd=0.5),cl=rep(cl1,nsim))
d2 <- data.frame(x=rnorm(n=nsim,sd=2),y=rnorm(n=nsim,sd=2),cl=rep(cl2,nsim))
d3 <- data.frame(x=c(rnorm(n=nsim/2,mean=-1,sd=0.5),rnorm(n=nsim/2,mean=1,sd=0.5)), y=  c(rnorm(n=nsim/2,mean=-1,sd=0.5),rnorm(n=nsim/2,mean=1,sd=0.5)), cl=rep(cl3,nsim))

pdf("/users/ndyjack/Dist_Proj/figures_final/figure_1.4E.pdf",width=6,height=4)
  plot.new()

  #Figure
  par(new = "TRUE",plt = c(0.05,0.30,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
  box()
  points(x=d1$x,y=d1$y,col=d1$cl,cex=0.8,pch=16)

  #Figure
  par(new = "TRUE",plt = c(0.35,0.60,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
  box()
  points(x=d2$x,y=d2$y,col=d2$cl,cex=0.8,pch=16)

  #Figure
  par(new = "TRUE",plt = c(0.65,0.90,0.45,0.80),las = 1,cex.axis = 1)
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
  box()
  points(x=d3$x,y=d3$y,col=d3$cl,cex=0.8,pch=16)

dev.off()


