#supplemental figure 1, Gap(k) for real data x each geneset
#scaled and unscaled for each
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)

map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

gsets = c(
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
)
gsets <- rev(gsets)
gsets <- sub("_","",gsets)

dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

k_tmp <- 1:5

results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/'
res_u <- lapply(gsets, function(x) sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,x,'.',y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  return(tmp)
}))
res_s <- lapply(res_u, function(y) apply(y,2, function(x) map_fxn(x,0,1,min(x),max(x))))

xlim <- c(0.9,5.1)
xticks <- k_tmp

ylim_s <- c(0.00,1.0)
yticks_s <- c(0.01,0.99) #c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_s <- c('min','max') #formatC(yticks_s,format='f',digits=1)

ylim_u <- c(0.00,3.0)
yticks_u <- c(0.00,1.00,2.00,3.00)
ylabs_u <- formatC(yticks_u,format='f',digits=1)

palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))
pch_labs <- rep(c(15,16,17),5)
pt_size <- 1.5


plot_locs <- list(
  list(c(0.05,0.28,0.50,0.85),c(0.05,0.28,0.07,0.42)),
  list(c(0.29,0.51,0.50,0.85),c(0.29,0.51,0.07,0.42)),
  list(c(0.52,0.74,0.50,0.85),c(0.52,0.74,0.07,0.42)),
  list(c(0.75,0.97,0.50,0.85),c(0.75,0.97,0.07,0.42))
)

pdf("/users/ndyjack/Dist_Proj/figures_final/supfig_1.1.pdf",width=12,height=9)

  #dist legend (top)
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')

for(j in 1:length(gsets)){
  #scaled lineplot
  par(new = "TRUE",plt = plot_locs[[j]][[1]])
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_s,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    points(x=k_tmp,y=res_s[[j]][,i],type='b',pch=pch_labs[i],col=col_labs[i],cex=pt_size)
  }
  mtext(side=1,text=gsets[j],cex=1.2,line=1.0,font=2,las=1)
  if(j==1){
    axis(side=2,labels=ylabs_s,at=yticks_s,cex=1.1,las=2,mgp=c(3, .5, 0))
    mtext(side=2,text="Gap Statistic (scaled)",cex=1.2,line=2.0,font=2,las=3)
  }

  #uncaled lineplot
  par(new = "TRUE",plt = plot_locs[[j]][[2]])
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim_u,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    points(x=k_tmp,y=res_u[[j]][,i],type='b',pch=pch_labs[i],col=col_labs[i],cex=pt_size)
  }
  axis(side=1,labels=xticks,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=1,text="k(# clusters)",cex=1.2,line=1.5,font=2,las=1)
  if(j==1){
    axis(side=2,labels=ylabs_u,at=yticks_u,cex=1.1,las=2,mgp=c(3, .5, 0))
    mtext(side=2,text="Gap Statistic (unscaled)",cex=1.2,line=2.0,font=2,las=3)
  }
}
dev.off()


res_wu <- do.call(rbind,res_u)
res_ws <- do.call(rbind,res_s) 
res_nmz <- as.vector(sapply(gsets, function(i) sapply(k_tmp, function(j) paste0(i,'_k=',j))))
rownames(res_ws) <- rownames(res_wu) <- res_nmz
write.csv(res_wu,"/users/ndyjack/Dist_Proj/tables_final/suptab_1.1.csv",quote=F)
write.csv(res_ws,"/users/ndyjack/Dist_Proj/tables_final/suptab_2.1.csv",quote=F)

