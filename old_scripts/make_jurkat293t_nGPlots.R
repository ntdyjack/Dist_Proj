#read in simulation results
in_dir <- "/users/ndyjack/Dist_Proj/rdata/"
sim_files <- list.files(in_dir)
sim_files <- sim_files[grep("jurkat293t.nG.[0-9]",sim_files)]
sim_files <- sim_files[order(as.numeric(sapply(sim_files, function(x) strsplit(x,"\\.")[[1]][3])))]
sim_files <- paste0(in_dir, sim_files)
sim_results <- lapply(sim_files, function(x) readRDS(x))

#plotting parameters
nG_vals <- unlist(lapply(sim_results, function(x) lapply(x, function(y) y$nGenes)))
xmin <- min(nG_vals)
xmax <- max(nG_vals) 
cols <- c("#0061FF","#2C9417","#F700FF")
fills <- paste0(cols,"35")
xticks <- pretty(c(xmin,xmax),5)
xlabs <- formatC(xticks,format="f",digits=0)
xaxis <- "nGenes (decreasing by variance)"
title_lab <- 'Jurkat(n=1607) and 293t(n=1781); 3799-55987 nUMIs (med=15041)'

qtls <- c(0.25,0.50,0.75)
plotdat_sil <- lapply(c('euc','poi','lrm'), function(x) {
  do.call(rbind,lapply(sim_results, function(y) t(sapply(y, function(z) quantile(z$sil[,x],qtls)))))
})
plotdat_ari <- lapply(c('euc','poi','lrm'), function(x) {
  unlist(lapply(sim_results, function(y) sapply(y, function(z) z$ari[x])))
})
plotdat_time <- lapply(c(1,3,2), function(x) {
  unlist(lapply(sim_results, function(y) sapply(y, function(z) z$time[x])))
})

tick_min <- 0.05
ymin_sil <- (round(min(unlist(plotdat_sil)),2)/tick_min)*tick_min
ymax_sil <- (round(max(unlist(plotdat_sil)),2)/tick_min)*tick_min
yticks_sil <- pretty(c(ymin_sil,ymax_sil),4)
ylabs_sil <- formatC(yticks_sil,format="f",digits=1)

ymin_ari <- (round(min(unlist(plotdat_ari)),2)/tick_min)*tick_min
ymax_ari <- (round(max(unlist(plotdat_ari)),2)/tick_min)*tick_min
yticks_ari <- pretty(c(ymin_ari,ymax_ari),4)
ylabs_ari <- formatC(yticks_ari,format="f",digits=1)

tick_min <- 100
ymin_time <- (round(min(unlist(plotdat_time)),0)/tick_min)*tick_min
ymax_time <- (round(max(unlist(plotdat_time)),0)/tick_min)*tick_min
yticks_time <- pretty(c(ymin_time,ymax_time),4)
ylabs_time <- formatC(yticks_time,format="f",digits=0)

yaxis_sil <- "Silhouette Score"
yaxis_ari <- "Adjusted Rand Index"
yaxis_time <- "time (seconds)"
leglabs <- c("Euclidean","Poisson","Multinomial") 

pdf("/users/ndyjack/Dist_Proj/images/jurkat293t_nG_tests.pdf",width=8,height=6)
par(mfrow=c(1,1),mar=c(3.0,4.0,2.0,1.0))
plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_sil),max(yticks_sil)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_sil,at=yticks_sil,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=2,font=2)
mtext(side=2,text=yaxis_sil,cex=1.0,line=3,font=2)
mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
legend("topright",legend=leglabs,col=cols,lty=1,lwd=1,cex=1,bty='n',pch=16)
for(i in 1:3){
  lines(x=nG_vals,y=plotdat_sil[[i]][,2],lty=1,col=cols[i],lwd=1)
  points(x=nG_vals,y=plotdat_sil[[i]][,2],pch=16,col=cols[i],cex=1)
  polygon(x=c(nG_vals,rev(nG_vals)), 
    y=c(plotdat_sil[[i]][,1],rev(plotdat_sil[[i]][,3])),border=NA,col=fills[i])
}

plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_ari),max(yticks_ari)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=2,font=2)
mtext(side=2,text=yaxis_ari,cex=1.0,line=3,font=2)
mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
legend("right",legend=leglabs,col=cols,lty=1,lwd=1,cex=1,bty='n',pch=16)
for(i in 1:3){
  lines(x=nG_vals,y=plotdat_ari[[i]],lty=1,col=cols[i],lwd=1)
  points(x=nG_vals,y=plotdat_ari[[i]],pch=16,col=cols[i],cex=1)
}

plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(min(yticks_time),max(yticks_time)),
  main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs_time,at=yticks_time,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=2,font=2)
mtext(side=2,text=yaxis_time,cex=1.0,line=3,font=2)
mtext(side=3,text=title_lab,cex=1.0,line=0.5,font=2)
legend("topleft",legend=leglabs,col=cols,lty=1,lwd=1,cex=1,bty='n',pch=16)
for(i in 1:3){
  lines(x=nG_vals,y=plotdat_time[[i]],lty=1,col=cols[i],lwd=1)
  points(x=nG_vals,y=plotdat_time[[i]],pch=16,col=cols[i],cex=1)
}
dev.off()


