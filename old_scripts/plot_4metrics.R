#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)
rndup <- function(x,d) { ceiling(x*d) / d }
rnddn <- function(x,d) { floor(x*d) / d }

#####
files = c(
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl",
  "cellbench10x_rnaseq_3cl",
  "cellbench10x_rnaseq_5cl",
  "cellbenchDS_rnaseq_3cl",
  "cellbenchCS2_rnaseq_3cl",
  "cellbenchCS2_rnaseq_5cl"
)

gsets = c(
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
)
gsets <- rev(gsets)

files_clvec <- as.vector(sapply(1:length(files), function(i) rep(i,4)))

datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5, 
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

#distance groupings
dists_grp1 <- c(
  'L1'='Geometric', 'L2'='Geometric', 'SQL2'='Geometric',
  'JSM'='Metric', 'JSD'='PseudoMetric', 'MKL'='BregmanDiv',
  'HEL'='Metric', 'BCM'='Metric', 'BCD'='Dissimilarity',
  'TVD'='Geometric', 'LLR'='Dissimilarity', 'UWLLR'='Dissimilarity',
  'ISD'='Probabilistic', 'RISD'='Probabilistic', 'SIS'='Probabilistic'
)

dists_grp2 <- c(
  'L1'='Geometric', 'L2'='Geometric', 'SQL2'='Geometric',
  'JSM'='Probabilistic', 'JSD'='Probabilistic', 'MKL'='Probabilistic',
  'HEL'='Probabilistic', 'BCM'='Probabilistic', 'BCD'='Probabilistic',
  'TVD'='Probabilistic', 'LLR'='Probabilistic', 'UWLLR'='Probabilistic',
  'ISD'='Probabilistic', 'RISD'='Probabilistic', 'SIS'='Probabilistic'
)



#  'L1'='Geometric', 'L2'='Geometric', 'SQRL2'='Geometric', 'JSM/PSM'='Probabilistic',
#  'JSD/PSD'='Probabilistic', 'MKL'='Probabilistic', #'POISSON'=6,
#  'HELLINGER'='Probabilistic',
#  'BHAT_METR'='Probabilistic', 'BHAT_DIST'='Probabilistic', 'TVD'='Probabilistic','LLR'='Probabilistic',
#  #'REV_MKL'=12, 'REV_POIS'=13 ,
#  'ITA_SAI'='Probabilistic',
  #'REV_ITA_SAI'=17,
#  'COS'='Geometric',
  #'COS_SIM'=20,'PRB_COS_SIM'=21,
#  'PL2'='Probabilistic','PSL2'='Probabilistic'#,
#  'std_pca'='std_pca','trn_pca'='trn_pca','glm_pca'='glm_pca'
#)

#load in results files,ARI
results_dir <- "/users/ndyjack/Dist_Proj/tables/ari_new/"
ari_vals <- results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".kmeds.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
ari_means <- dist_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))

#k accuracy
results_dir <- "/users/ndyjack/Dist_Proj/tables/kacc_new/"
kacc_vals <- results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".k30.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
kacc_means <- dist_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))


#G+
results_dir <- "/users/ndyjack/Dist_Proj/tables/gplus_new/"
gplus_sorted <- results_sorted <- 1 - sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
gplus_means <- dist_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))

#trajcor
files = c(
  "cellmix_rnaseq_traj",
  "rnamixcs2_rnaseq_traj",
  "rnamixss_rnaseq_traj"
)
datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

results_dir <- "/users/ndyjack/Dist_Proj/tables/trajcor_new/"
cor_sorted <- results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- read.table(paste0(results_dir,y,".",x,".txt"),header=T)
  return(mean(tmp$pearson))
}))
cor_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))


metric_nmz <- c("ARI","kAcc","1-G+","RankCor")
#metric_names <- c("ARI","kAcc","G+","Cor")

metrics_by_gset <- lapply(gsets, function(y) {
  x <- sub("_","",y)
  tmp <- cbind(ari_means[x,],kacc_means[x,],gplus_means[x,],cor_means[x,])
  colnames(tmp) <- metric_nmz
  rownames(tmp) <- names(dists_vec)
  return(tmp)
})

metrics_avgs <- as.data.frame(cbind(colMeans(ari_means),colMeans(kacc_means),colMeans(gplus_means),colMeans(cor_means)))
colnames(metrics_avgs) <- metric_nmz
rownames(metrics_avgs) <- names(dists_vec) 

metrics_by_gset[[5]] <- metrics_avgs
gsets[5] <- "_avg"

palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
cols <- sample(palette3_all, length(dists_vec))

#plotting parameters
zlim <- c(0,1)
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
nrow <- length(dists_vec)
ncol <- 4
xlim <- c(0,ncol+2)
ylim <-c(0,nrow+1)
zleglocs <- c(0.01,0.25,0.50,0.75,0.99)
zleglabs <- c("0.00","0.25","0.50","0.75","1.00")
#metric_nmz <- c("ARI","kAcc","1-G+","RankCor")
dist_gcol <- ifelse(dists_grp=="Probabilistic",'red','blue')
dist_pch <- ifelse(dists_grp=="Probabilistic",15,16)

pdf("/users/ndyjack/Dist_Proj/images/metrics_heatmap_unscaled.pdf",width=8,height=12)
for(i in 5:1){
  #draw heatmap
  tmp <- metrics_by_gset[[i]]
  o <- order(rowMeans(tmp),decreasing=F)
  tmp <- tmp[o,]
  plot.new()
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.80),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA) 
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=nrow+0.7,labels=metric_nmz)
  text(x=ncol+1.0,y=(1:nrow),labels=rownames(tmp),col=dist_gcol[o])

  #color legend
  par(new = "TRUE",plt = c(0.15,0.65,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))

  #dist type legend
  par(new = "TRUE",plt = c(0.70,0.90,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=rep(0.5,2),y=c(0.2,0.8),labels=c("Probabilistic","Geometric"),col=c("red","blue"),cex=1.5)

  #title
  par(new = "TRUE",plt = c(0.30,0.70,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.5)

}
dev.off()

zleglocs <- c(0.01,0.99)
zleglabs <- c("min","max")

pdf("/users/ndyjack/Dist_Proj/images/metrics_heatmap_scaled.pdf",width=8,height=12)
for(i in 5:1){
  #draw heatmap
  tmp <- metrics_by_gset[[i]]
  o <- order(rowMeans(tmp),decreasing=F)
  tmp <- tmp[o,]
  plot.new()
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.80),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:ncol){
    zlim_tmp <- c(min(tmp[,k]),max(tmp[,k]))
    col_pal_tmp <- colorNumeric(col_pal, domain=zlim_tmp, na.color = "#808080",alpha = FALSE,reverse = FALSE)
    for(j in 1:nrow){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_tmp(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=nrow+0.7,labels=metric_nmz)
  text(x=ncol+1.0,y=(1:nrow),labels=rownames(tmp),col=dist_gcol[o])

  #color legend
  par(new = "TRUE",plt = c(0.15,0.65,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))

  #dist type legend
  par(new = "TRUE",plt = c(0.70,0.90,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=rep(0.5,2),y=c(0.2,0.8),labels=c("Probabilistic","Geometric"),col=c("red","blue"),cex=1.5)

  #title
  par(new = "TRUE",plt = c(0.30,0.70,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.5)
}
dev.off()




zlim_ari <- c(
  min(sapply(metrics_by_gset, function(x) min(x[,'ARI']))),
  max(sapply(metrics_by_gset, function(x) max(x[,'ARI'])))
)
col_pal_ari <- colorNumeric(col_pal, domain=zlim_ari, na.color = "#808080",alpha = FALSE,reverse = FALSE)

zlim_kacc <- c(
  min(sapply(metrics_by_gset, function(x) min(x[,'kAcc']))),
  max(sapply(metrics_by_gset, function(x) max(x[,'kAcc'])))
)
col_pal_kacc <- colorNumeric(col_pal, domain=zlim_kacc, na.color = "#808080",alpha = FALSE,reverse = FALSE)

zlim_gplus <- c(
  min(sapply(metrics_by_gset, function(x) min(x[,'G+']))),
  max(sapply(metrics_by_gset, function(x) max(x[,'G+'])))
)
col_pal_gplus <- colorNumeric(col_pal, domain=zlim_gplus, na.color = "#808080",alpha = FALSE,reverse = FALSE)

zlim_cor <- c(
  min(sapply(metrics_by_gset, function(x) min(x[,'Cor']))),
  max(sapply(metrics_by_gset, function(x) max(x[,'Cor'])))
)
col_pal_cor <- colorNumeric(col_pal, domain=zlim_cor, na.color = "#808080",alpha = FALSE,reverse = FALSE)

#zleglocs <- c(0.01,0.99)
#zleglabs <- c("min","max")

xlim <- c(0.1,0.9) #ari
xloc <- pretty(zlim_ari,4)
xlab <- formatC(xloc,digits=1,format='f')
ylim <- c(0.6,1.0)
yloc <- pretty(ylim,4)
ylab <- formatC(yloc,digits=1,format='f')

cex_min <- 1.0
cex_max <- 4.0
#output = output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start)
map_cor <- function(x) {
  cex_min + ((cex_max - cex_min) / (zlim_cor[2] - zlim_cor[1]) * (x - zlim_cor[1] ) )
}
trns_min <- 100
trns_max <- 255
map_trns <- function(x) {
  trns_min + ((trns_max - trns_min) / (zlim_gplus[2] - zlim_gplus[1]) * (x - zlim_gplus[1] ) )
}

pdf("/users/ndyjack/Dist_Proj/images/metrics_plots_scatterplot.pdf",width=8,height=12)
for(i in 4:1){
  #dist legend
  plot.new()
  par(new = "TRUE",plt = c(0.05,0.95,0.65,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('topleft',legend=names(dists_vec),col=cols,ncol=5,
    x.intersp=0.5,y.intersp=2,text.width=0.17,cex=1.2,
    pch=16,pt.cex=2.0,bty='n')

  #title
  par(new = "TRUE",plt = c(0.30,0.70,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.5)

  #transparcency and size legends
  par(new = "TRUE",plt = c(0.05,0.95,0.65,0.75),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=c(0.20,0.35,0.65,0.80),y=rep(0.5,4),labels=c("low Cor","high Cor","low G+","high G+"),cex=1.5)
  points(x=c(0.20,0.35,0.65,0.80),y=rep(0.7,4),cex=c(1.0,4.0,2.5,2.5),
    pch=16,col=c('#000000','#000000','#00000064','#000000FF'))

  #scatter plot
  tmp <- metrics_by_gset[[i]]
  xval <- tmp[,'ARI']
  yval <- tmp[,'kAcc']
  cval <- tmp[,'Cor']
  tval <- tmp[,'G+']
  cex_vals <- map_cor(cval)
  trns_vals <- as.hexmode(round(map_trns(tval)))
  col_vals <- paste0(cols,trns_vals)
  par(new = "TRUE",plt = c(0.15,0.95,0.10,0.65),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
  box()
  points(x=xval,y=yval,col=col_vals,pch=16,cex=cex_vals)
  axis(side=1,labels=xlab,at=xloc,cex=0.7,las=1,mgp=c(1.0, .4, 0))
  mtext(side=1,text='ARI',cex=1.5,line=1.5,las=1)
  mtext(side=2,text='kAcc',cex=1.5,line=2.5,las=2)
  axis(side=2,labels=ylab,at=yloc,cex=0.7,las=2,mgp=c(3, .7, 0))
}
dev.off()

#all pairs of metrics
metrics_pairs <- combn(metric_nmz,2)
#c("ARI","knnAcc","1-G+","RankCor") 

sixplot_lims <- list(
  c(0.07,0.31,0.48,0.82),
  c(0.40,0.64,0.48,0.82),
  c(0.73,0.97,0.48,0.82),
  c(0.07,0.31,0.07,0.41),
  c(0.40,0.64,0.07,0.41),
  c(0.73,0.97,0.07,0.41)
)


pdf("/users/ndyjack/Dist_Proj/images/metrics_2x3_scatterplots.pdf",width=12,height=8)
for(i in 5:1){
  #dist legend
  plot.new()
  par(new = "TRUE",plt = c(0.05,1.00,0.75,0.94),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('topleft',legend=names(dists_vec),col=cols,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=dist_pch,pt.cex=2.0,bty='n')
  points(x=rep(0.85,2),y=c(0.60,0.80),col='black',cex=2,pch=c(15,16))
  text(x=rep(0.85,2),y=c(0.60,0.80),labels=c("Prob","Geom"),cex=1.2,pos=4)

  #title
  par(new = "TRUE",plt = c(0.30,0.70,0.92,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.5)

  #3x2 scatterplot
  for(j in 1:ncol(metrics_pairs)){
    par(new='TRUE',las=1,cex.axis=1,plt = sixplot_lims[[j]],bg="transparent")
    #plot.new() 
    xlab <- metrics_pairs[1,j]
    ylab <- metrics_pairs[2,j]
    tmp <- metrics_by_gset[[i]]
    xval <- tmp[,xlab]
    yval <- tmp[,ylab]
    xmin <- rnddn(min(xval),10)
    xmax <- rndup(max(xval),10)
    xlim <- c(xmin,xmax)
    xloc <- pretty(xlim,4)
    xtxt <- formatC(xloc,digits=2,format='f')
    ymin <- rnddn(min(yval),10)
    ymax <- rndup(max(yval),10)
    ylim <- c(ymin,ymax)
    yloc <- pretty(ylim,4)
    ytxt <- formatC(yloc,digits=2,format='f')
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
    box()
    points(x=xval,y=yval,col=cols,cex=1.5,pch=dist_pch)
    axis(side=1,labels=xtxt,at=xloc,cex=0.7,las=1,mgp=c(1.0, .4, 0))
    mtext(side=1,text=xlab,cex=1.2,line=1.5,las=1)
    mtext(side=2,text=ylab,cex=1.2,line=2.5,las=3)
    axis(side=2,labels=ytxt,at=yloc,cex=0.7,las=2,mgp=c(3, .7, 0))
  }     
}
dev.off()


