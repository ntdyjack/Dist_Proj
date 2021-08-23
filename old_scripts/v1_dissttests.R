#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman")}
library(pacman)
p_load(Matrix,cluster,mclust,Rfgc,NbClust)
#library(Matrix)
#library(cluster)
#library(mclust)
#library(Rfgc)
#library(NbClust)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
k_val <- as.integer(sub("cl","",strsplit(args[3],split="_")[[1]][3]))

#location of data file
#cluster_lab <-  as.numeric(readLines("/users/ndyjack/Dist_Proj/tables/jurkat293t_2cl_labs.txt"))
#counts_loc <- "/users/ndyjack/Dist_Proj/tables/jurkat293t_2cl_expr.mtx"
counts_loc <- paste0(args[1],args[3],"_expr.mtx") 
counts_mat <- readMM(counts_loc)
counts_mat <- t(as.matrix(counts_mat))
cluster_lab <- as.numeric(readLines(paste0(args[1],args[3],"_labs.txt")))
counts_mat <- counts_mat[1:500,]
cluster_lab <- cluster_lab[1:500]

#names of the distances being tested
dists_vec <- c('L1','L2','SQRL2','JSM/PSM','JSD/PSD','MKL','POISSON','HELLINGER',
  'BHAT_METRIC','BHAT_DISTANCE','TVD','LLR',
  'EMD','REV_MKL','REV_POISSON','UWLLR','OLLR')
#idx_test <- 1:12 
idx_test <- 1:length(dists_vec)

nCell_vec <- floor(seq(.05,1.0,length.out=50)*length(cluster_lab))
qtls <- c(0.40,0.50,0.60)

results_list <- lapply(idx_test, function(i) {
  cat("running... ",dists_vec[i],"\n")
  start_time <- Sys.time()
#  dist_tmp <- as.matrix(dist(counts_mat))
  dist_tmp <- Rfgc::dist_matrixdf(counts_mat, i-1)
  time_tmp <- Sys.time() - start_time

  stats_tmp <- lapply(1:length(nCell_vec), function(j) {
    set.seed(1752)
    cells_use <- sample(1:length(cluster_lab),size=nCell_vec[j])
    labs_use <- cluster_lab[cells_use]
    dist_tmp_use <- as.dist(dist_tmp[cells_use,cells_use])
    sil_tmp <- silhouette(x=labs_use,dist=dist_tmp_use)[,3]
    sil_tmp <- quantile(sil_tmp, qtls)
#    sil_tmp <- c(mean(sil_tmp) - sd(sil_tmp),mean(sil_tmp), mean(sil_tmp) + sd(sil_tmp)) 
    ari_tmp <- adjustedRandIndex(pam(dist_tmp_use,k=k_val)$cluster,labs_use)
    return(list(sil=sil_tmp,ari=ari_tmp))
  })
  
  return(list(
    sil = sapply(stats_tmp, function(x) x$sil),
    ari = sapply(stats_tmp, function(x) x$ari),
    time = time_tmp
  ))

})

cat("plotting results...\n")

#plotting parameters
xmin <- min(nCell_vec)
xmax <- max(nCell_vec)
cols <- gg_color_hue(length(results_list))
fills <- paste0(cols,"35")
shps <- rep(21:25,5)[idx_test]
ltys <- rep(1:5,5)[idx_test]
xticks <- pretty(c(xmin,xmax),5)
xlabs <- formatC(xticks,format="f",digits=0)
xaxis <- "nCell"
title_lab <- args[3]

tick_min <- 0.05
ymin_sil <- (round(min(sapply(results_list, function(x) x$sil)),2)/tick_min)*tick_min
ymax_sil <- (round(max(sapply(results_list, function(x) x$sil)),2)/tick_min)*tick_min
yticks_sil <- pretty(c(ymin_sil,ymax_sil),4)
ylabs_sil <- formatC(yticks_sil,format="f",digits=1)

ymin_ari <- (round(min(sapply(results_list, function(x) x$ari)),2)/tick_min)*tick_min
ymax_ari <- (round(max(sapply(results_list, function(x) x$ari)),2)/tick_min)*tick_min
yticks_ari <- pretty(c(ymin_ari,ymax_ari),4)
ylabs_ari <- formatC(yticks_ari,format="f",digits=1)

tick_min <- 1
ymin_time <- (round(min(sapply(results_list, function(x) x$time)),0)/tick_min)*tick_min
ymax_time <- (round(max(sapply(results_list, function(x) x$time)),0)/tick_min)*tick_min
yticks_time <- pretty(c(ymin_time,ymax_time),4)
ylabs_time <- formatC(yticks_time,format="f",digits=0)

yaxis_sil <- "Silhouette Score"
yaxis_ari <- "Adjusted Rand Index"
yaxis_time <- "Time (seconds)"

draw_legend <- function(points=T){
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.23),las = 1, cex.axis = 1)
  plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n", yaxt="n",bty='n')
  if(points){
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=3,
      x.intersp=0.6,y.intersp=1,text.width=0.2,cex=1.0,
      pch=shps,pt.bg="#00000050",pt.cex=1.0,bty='n')
  }else{
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=3,
      x.intersp=0.5,y.intersp=1,text.width=0.2,cex=1.0,
      lty=ltys,lwd=1.5,bty='n')  
  }
}


pdf(paste0(args[2],args[3],"_disttest.pdf"),width=8,height=6)

plot.new()
par(new = "TRUE",plt = c(0.10,0.95,0.30,0.95))
plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_sil),max(yticks_sil)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_sil,at=yticks_sil,las=2,cex=0.8,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=1.5,font=2)
mtext(side=2,text=yaxis_sil,cex=1.0,line=2.0,font=2,las=3)
mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
for(i in 1:length(results_list)){
  lines(x=nCell_vec,y=results_list[[i]]$sil[2,],lty=ltys[i],col=cols[i],lwd=1.5)
  polygon(x=c(nCell_vec,rev(nCell_vec)),
    y=c(results_list[[i]]$sil[1,],rev(results_list[[i]]$sil[3,])),border=NA,col=fills[i])
}
draw_legend(F)

plot.new()
par(new = "TRUE",plt = c(0.10,0.95,0.30,0.95))
plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_ari),max(yticks_ari)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_ari,at=yticks_ari,las=2,cex=0.8,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=1.5,font=2)
mtext(side=2,text=yaxis_ari,cex=1.0,line=2.0,font=2,las=3)
mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
for(i in 1:length(results_list)){
  lines(x=nCell_vec,y=results_list[[i]]$ari,lty=ltys[i],col=cols[i],lwd=1.5)
}
draw_legend(F)

plot.new()
par(new = "TRUE",plt = c(0.10,0.95,0.30,0.95))
plot(x=1,y=1,type='n',
  xlim=c(min(yticks_time),max(yticks_time)),
  ylim=c(min(yticks_sil),max(yticks_sil)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n")
axis(side=1,labels=ylabs_time,at=yticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_sil,at=yticks_sil,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=yaxis_time,cex=1.0,line=1.5,font=2)
mtext(side=2,text=yaxis_sil,cex=1.0,line=2.0,font=2,las=3)
mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
for(i in 1:length(results_list)){
  points(x=results_list[[i]]$time,y=tail(results_list[[i]]$sil[2,],1),
    pch=shps[i],col=cols[i],bg="#00000050",cex=2)
}
draw_legend(T)


plot.new()
par(new = "TRUE",plt = c(0.10,0.95,0.30,0.95))
plot(x=1,y=1,type='n',
  xlim=c(min(yticks_time),max(yticks_time)),
  ylim=c(min(yticks_ari),max(yticks_ari)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n")
axis(side=1,labels=ylabs_time,at=yticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=yaxis_time,cex=1.0,line=1.5,font=2)
mtext(side=2,text=yaxis_ari,cex=1.0,line=2.0,font=2,las=3)
mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
for(i in 1:length(results_list)){
  points(x=results_list[[i]]$time,y=tail(results_list[[i]]$ari[2,],1),
    pch=shps[i],col=cols[i],bg="#00000050",cex=2)
}
draw_legend(T)

dev.off()


#plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_sil),max(yticks_sil)),
#  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
#axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
#axis(side=2,labels=ylabs_sil,at=yticks_sil,cex=0.8,las=2,mgp=c(3, .5, 0))
#mtext(side=1,text=xaxis,cex=1.0,line=2,font=2)
#mtext(side=2,text=yaxis_sil,cex=1.0,line=3,font=2)
#mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
#legend("topleft",legend=dists_vec[idx_test],col=cols,lty=1,lwd=2,bty='n',
#   x.intersp=0.5,text.width=nCell_vec[1],ncol=length(results_list)/2,cex=0.5)
#for(i in 1:length(results_list)){
#  lines(x=nCell_vec,y=results_list[[i]]$sil[2,],lty=1,col=cols[i],lwd=2)
  #points(x=nCell_vec,y=plotdat_sil[[i]][2,],pch=16,col=cols[i],cex=1)
#  polygon(x=c(nCell_vec,rev(nCell_vec)),
#    y=c(results_list[[i]]$sil[1,],rev(results_list[[i]]$sil[3,])),border=NA,col=fills[i])
#}

#plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_ari),max(yticks_ari)),
#  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
#axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
#axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=0.8,las=2,mgp=c(3, .5, 0))
#mtext(side=1,text=xaxis,cex=1.0,line=2,font=2)
#mtext(side=2,text=yaxis_ari,cex=1.0,line=3,font=2)
#mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
#legend("topleft",legend=dists_vec[idx_test],col=cols,lty=1,lwd=2,bty='n',
#  ncol=length(results_list)/2,x.intersp=0.5,text.width=nCell_vec[1],cex=0.5)
#for(i in 1:length(results_list)){
#  lines(x=nCell_vec,y=results_list[[i]]$ari,lty=1,col=cols[i],lwd=2)
  #points(x=nCell_vec,y=plotdat_ari[[i]],pch=16,col=cols[i],cex=1)
#}

#plot(x=1,y=1,type='n',
#  xlim=c(min(yticks_time),max(yticks_time)),
#  ylim=c(min(yticks_ari),max(yticks_ari)),
#  main="",xlab="",ylab="",xaxt="n", yaxt="n")
#axis(side=1,labels=ylabs_time,at=yticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
#axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=0.8,las=2,mgp=c(3, .5, 0))
#mtext(side=1,text=yaxis_time,cex=1.0,line=2,font=2)
#mtext(side=2,text=yaxis_ari,cex=1.0,line=3,font=2)
#mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
#legend("topleft",legend=dists_vec[idx_test],col=cols,bty='n',
#  ncol=length(results_list)/2,x.intersp=0.6,text.width=ymin_time/tick_min,cex=0.5,pch=16,pt.cex=1)
#for(i in 1:length(results_list)){
  #lines(x=nCell_vec,y=results_list[[i]]$ari,lty=1,col=cols[i],lwd=2)
#  points(x=results_list[[i]]$time,y=tail(results_list[[i]]$ari,1),
#    pch=16,col=cols[i],cex=2)
#}

#plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n", yaxt="n",bty='n')
#legend('topleft',legend=dists_vec[idx_test],col=cols,ncol=2,
#  x.intersp=0.6,y.intersp=1,text.width=0.3,cex=1.5,
#  pch=16,pt.cex=2,bty='n')
#mtext(side=3,text='Color Legend',cex=1.0,line=0.5,font=2)
