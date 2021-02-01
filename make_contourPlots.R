#read in simulation results
in_dir <- "/users/ndyjack/Dist_Proj/rdata/"
sim_files <- list.files(in_dir)
sim_files <- sim_files[grep("splatterTests.lF",sim_files)]
#sim_files <- sim_files[c(6,5,8,7,10,9,2,1,4,3)]
sim_files <- paste0(in_dir, sim_files)
#sim_files <- sim_files[order(as.numeric(sapply(sim_files, function(x) strsplit(x,"\\.")[[1]][[3]])))]
sim_results <- lapply(sim_files, function(x) readRDS(x))



test_info <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$info))))
test_info <- do.call(rbind,test_info)




euc_results <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$sil$euc))))
euc_results <- do.call(rbind,euc_results)

lrm_results <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$sil$lrm))))
lrm_results <- do.call(rbind,lrm_results)

poi_results <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$sil$poi))))
poi_results <- do.call(rbind,poi_results)

#plotting params
nticks <- 5
zmin <- 0 #-.1 #signif(min(euc_results,lrm_results,poi_results),1) - 0.1
zmax <- 1
zlim <- c(zmin,zmax)
ztickloc <- seq(zmin,zmax,length.out=nticks)
zticklab <- formatC(ztickloc,format="f",digits=1)
nlevs <- 100
levs <- pretty(zlim,nlevs)
library(viridis)
col_pal <- viridis(length(levs)-1)

#plotting coordinates (just the indices of the data matrices)
xcords <- ycords <- 1:nrow(euc_results)
xlim <- ylim <- c(min(xcords),max(xcords))
#xcords <- 1:nrow(euc_results)
#xlim <- c(min(xcords),max(xcords))
#ycords <- 1:ncol(euc_results)
#ylim <- c(min(ycords),max(ycords))

#axis ticks/ text labels
xtickloc <- ytickloc <- c(1,seq(5,20,length.out=nticks-1)) #ipretty(c,nticks)
#ytickloc <- c(1,seq(5,18,length.out=nticks-1))
xticks <- seq(0,9,length.out=nticks)
xticklab <- formatC(xticks,format="f",digits=1)
xaxis <- expression(bold("log" [2] * "fold change"))
yticks <- c(2.5,seq(25,100,length.out=nticks-1)) #seq(2.5,100,length.out=nticks)
yticklab <- formatC(yticks,format="f",digits=0)
yaxis <- "Percent Dropout"
xlab  <- "effect size for DEGs"
ylab <- "library size"

#want this to be a 3-panel plot with a legend (colorbar) on the far right
pdf("/users/ndyjack/Dist_Proj/images/splatterSim_Sil.pdf",width=12,height=4)
#initate new frame
plot.new()
#First plot (far left)
par(new = "TRUE",plt = c(0.08,0.30,0.17,0.88),las = 1, cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, euc_results, levels = levs, col = col_pal)
box()
mtext(side=3,text="Euclidean Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=-4,y=10,labels=yaxis,srt=90,cex=1.5,font=2)

#second plot (middle)
par(new = "TRUE",plt = c(0.36,0.58,0.17,0.88),las = 1,cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, poi_results, levels = levs, col = col_pal)
box()
mtext(side=3,text="Poisson Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=10,y=-2.5,labels=xaxis,font=2,cex=1.5)

#third plot (far right)
par(new = "TRUE",plt = c(0.64,0.86,0.17,0.88),las = 1,cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, lrm_results, levels = levs, col = col_pal)
box()
mtext(side=3,text="Multinomial Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)

#Add legend
par(new = "TRUE",plt = c(0.90,0.92,0.20,0.85),las = 1,cex.axis = 1)
plot.window(xlim = c(0, 1), ylim = range(levs), xaxs = "i",yaxs = "i")
rect(0, levs[-length(levs)], 1, levs[-1L], col = col_pal,border=NA)
box()
axis(side=4,labels=zticklab,at=ztickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=3.8,y=0.5,labels="Silhouette Score",cex=1.5,srt=270,font=2)
dev.off()


euc_times <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$times[1]))))
euc_times <- do.call(rbind,euc_times)

lrm_times <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$times[2]))))
lrm_times <- do.call(rbind,lrm_times)

poi_times <- lapply(sim_results, function(x) sapply(x, function(y) sapply(y, function(z) median(z$times[3]))))
poi_times <- do.call(rbind,poi_times)

zmin <- min(floor(cbind(euc_times,lrm_times,poi_times)))
zmax <- max(ceiling(cbind(euc_times,lrm_times,poi_times)))
zlim <- c(zmin,zmax)
ztickloc <- seq(zmin,zmax,length.out=nticks)
zticklab <- formatC(ztickloc,format="G",digits=1)
nlevs <- 100
levs <- pretty(zlim,nlevs)
col_pal <- viridis(length(levs)-1)

#want this to be a 3-panel plot with a legend (colorbar) on the far right
pdf("/users/ndyjack/Dist_Proj/images/splatterSim_Time.pdf",width=12,height=4)
#initate new frame
plot.new()
#First plot (far left)
par(new = "TRUE",plt = c(0.08,0.30,0.17,0.88),las = 1, cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, euc_times, levels = levs, col = col_pal)
box()
mtext(side=3,text="Euclidean Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=-4,y=10,labels=yaxis,srt=90,cex=1.5,font=2)

#second plot (middle)
par(new = "TRUE",plt = c(0.36,0.58,0.17,0.88),las = 1,cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, poi_times, levels = levs, col = col_pal)
box()
mtext(side=3,text="Poisson Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=10,y=-2.5,labels=xaxis,font=2,cex=1.5)

#third plot (far right)
par(new = "TRUE",plt = c(0.64,0.86,0.17,0.88),las = 1,cex.axis = 1)
plot.new()
plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
.filled.contour(xcords, ycords, lrm_times, levels = levs, col = col_pal)
box()
mtext(side=3,text="Multinomial Distance",font=2,line=0.1,cex=1.5)
axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)

#Add legend
par(new = "TRUE",plt = c(0.90,0.92,0.20,0.85),las = 1,cex.axis = 1)
plot.window(xlim = c(0, 1), ylim = range(levs), xaxs = "i",yaxs = "i")
rect(0, levs[-length(levs)], 1, levs[-1L], col = col_pal,border=NA)
box()
axis(side=4,labels=zticklab,at=ztickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
par(xpd=T)
text(x=4.0,y=median(levs),labels="Time (seconds)",cex=1.5,srt=270,font=2)
dev.off() 
