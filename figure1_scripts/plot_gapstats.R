args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)
rndup <- function(x,d) { ceiling(x*d) / d }
rnddn <- function(x,d) { floor(x*d) / d }


gsets = c(
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
)
gsets <- rev(gsets)
gsets <- sub("_","",gsets)

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

k_tmp <- 1:5

results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/'
results_list <- lapply(gsets, function(x) sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,x,'.',y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  return(tmp) 
}))

results_sims <- lapply(gsets, function(x) sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,'splat.',x,'.',y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  return(tmp)
}))

results_noise <- readLines(paste0(results_dir,'randdist.txt'))
results_noise <- as.numeric(strsplit(results_noise," ")[[1]])
results_noiss <- results_noise + abs(min(results_noise))
results_noiss <- results_noiss/max(results_noiss)
#results_noiss <- results_noise + (results_noise - min(results_noise))/(max(results_noise) - min(results_noise))
#results_noiss[2] <- ceiling(results_noiss[2])
#results_noiss[5] <- floor(results_noiss[5]) 
#results_noise/max(abs(results_noise))

xlim <- c(0.5,5.5) 
xticks <- k_tmp

ylim_s <- c(0.00,1.0)
yticks_s <- c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_s <- formatC(yticks_s,format='f',digits=1)

ylim_u <- c(-0.0001,3.0)
yticks_u <- c(0.00,1.00,2.00,3.00)
ylabs_u <- formatC(yticks_u,format='f',digits=1)


palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))
pch_labs <- rep(c(15,16,17),5)

pdf("/users/ndyjack/Dist_Proj/images/gapstat_k_real.pdf",width=12,height=12)
for(j in 1:length(results_list)){
  #dist legend (top)
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')


  #scaled lineplot
  par(new = "TRUE",plt = c(0.10,0.90,0.50,0.90))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_s,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    tmp <- results_list[[j]][,i]
    tmp <- tmp/max(tmp)
    points(x=k_tmp,y=tmp,type='b',pch=pch_labs[i],col=col_labs[i])
  }
  #points(x=k_tmp,y=results_noiss,type='b',pch=0,col='black')
  axis(side=2,labels=ylabs_s,at=ylabs_s,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticks,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text="Gap Statistic (scaled)",cex=1.2,line=2.5,font=2,las=3)
  mtext(side=3,text=gsets[j],cex=1.2,line=0.2,font=2,las=1)

  #uncaled lineplot
  par(new = "TRUE",plt = c(0.10,0.90,0.05,0.45))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_u,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    tmp <- results_list[[j]][,i]
    points(x=k_tmp,y=tmp,type='b',pch=pch_labs[i],col=col_labs[i])
  }
  #points(x=k_tmp,y=results_noise,type='b',pch=0,col='black')
  axis(side=2,labels=ylabs_u,at=yticks_u,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticks,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=1,text="k( # clusters)",cex=1.2,line=1.5,font=2,las=1)
  mtext(side=2,text="Gap Statistic (unscaled)",cex=1.2,line=2.5,font=2,las=3)
}
dev.off()


ylim_s <- c(0.90,1.0)
yticks_s <- c(0.90,0.95,1.0)
ylabs_s <- formatC(yticks_s,format='f',digits=1)

ylim_u <- c(-0.0001,2.0)
yticks_u <- c(0.00,1.00,2.00,3.00)
ylabs_u <- formatC(yticks_u,format='f',digits=1)


pdf("/users/ndyjack/Dist_Proj/images/gapstat_k_simulated.pdf",width=12,height=12)
for(j in 1:length(results_list)){
  #dist legend (top)
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')


  #scaled lineplot
  par(new = "TRUE",plt = c(0.10,0.90,0.50,0.90))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_s,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    tmp <- results_sims[[j]][,i]
    tmp <- tmp/max(tmp)
    points(x=k_tmp,y=tmp,type='b',pch=pch_labs[i],col=col_labs[i])
  }
  #points(x=k_tmp,y=results_noiss,type='b',pch=0,col='black')
  axis(side=2,labels=ylabs_s,at=ylabs_s,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticks,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text="Gap Statistic (scaled)",cex=1.2,line=2.5,font=2,las=3)
  mtext(side=3,text=gsets[j],cex=1.2,line=0.2,font=2,las=1)

  #uncaled lineplot
  par(new = "TRUE",plt = c(0.10,0.90,0.05,0.45))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_u,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    tmp <- results_sims[[j]][,i]
    points(x=k_tmp,y=tmp,type='b',pch=pch_labs[i],col=col_labs[i])
  }
  #points(x=k_tmp,y=results_noise,type='b',pch=0,col='black')
  axis(side=2,labels=ylabs_u,at=yticks_u,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticks,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=1,text="k( # clusters)",cex=1.2,line=1.5,font=2,las=1)
  mtext(side=2,text="Gap Statistic (unscaled)",cex=1.2,line=2.5,font=2,las=3)
}
dev.off()

