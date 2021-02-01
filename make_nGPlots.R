#read in simulation results
in_dir <- "/users/ndyjack/Dist_Proj/rdata/"
sim_files <- list.files(in_dir)
sim_files <- sim_files[grep("splatterTests.nG.[0-9]",sim_files)]
sim_files <- paste0(in_dir, sim_files)
sim_results <- lapply(sim_files, function(x) readRDS(x))

qtls <- c(0.25,0.50,0.75)

plotdat <- list(
  sapply(sim_results, function(x) sapply(x, function(y) y$times[1])),
  sapply(sim_results, function(x) sapply(x, function(y) y$times[3])),
  sapply(sim_results, function(x) sapply(x, function(y) y$times[2]))
)
plotdat <- lapply(plotdat, function(x) apply(x,1,function(y) quantile(y,qtls)))
#plotdat <- lapply(1:3, function(i) t(sapply(sim_results, function(x) quantile(sapply(x, function(y) y$time[i]),qtls) )))
#plotdat <- list(
#  lapply(sim_results, function(x) sapply(x, function(y) y$sil$euc)),
#  lapply(sim_results, function(x) sapply(x, function(y) y$sil$poi)),
#  lapply(sim_results, function(x) sapply(x, function(y) y$sil$lrm))
#)
#plotdat <- list(
#  sapply(1:9, function(i) quantile(sapply(plotdat[[1]], function(x) x[,i]),qtls)),
#  sapply(1:9, function(i) quantile(sapply(plotdat[[2]], function(x) x[,i]),qtls)),
#  sapply(1:9, function(i) quantile(sapply(plotdat[[3]], function(x) x[,i]),qtls))
#) 


#plotting parameters
xmin <- 100
xmax <- 20000
nG_vals <- c(100,500,1000,2000,3000,6000,10000,15000,20000)
#nG_vals <- c(100,500,1000,3000,5000,10000,20000)
cols <- c("#0061FF","#2C9417","#F700FF")
fills <- paste0(cols,"35")
xticks <- nG_vals
xlabs <- formatC(xticks,format="f",digits=0)
xaxis <- "nGenes"
ymin <- 0
ymax <- 6000
yticks <- pretty(c(ymin,ymax),5)
ylabs <- formatC(yticks,format="f",digits=1)
yaxis <- "Time (seconds)"
leglabs <- c("Euclidean","Poisson","Multinomial") 

pdf("/users/ndyjack/Dist_Proj/images/splatterSim_nG.pdf",width=12,height=6)
par(mfrow=c(1,1),mar=c(3.1,5.6,1.0,3.1))
plot(x=1,y=1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),main="",xlab="",ylab="",xaxt="n", yaxt="n",type='n')
axis(side=1,labels=xlabs,at=xticks,cex=0.8,las=1,mgp=c(3, .5, 0))
axis(side=2,labels=ylabs,at=yticks,cex=0.8,las=2,mgp=c(3, .5, 0))
mtext(side=1,text=xaxis,cex=1.0,line=2)
mtext(side=2,text=yaxis,cex=1.0,line=4)
legend("topleft",legend=leglabs,col=cols,lty=1:3,lwd=3,cex=1,bty='n')
for(i in 1:3){
	lines(x=nG_vals,y=plotdat[[i]][2,],lty=i,lwd=2,col=cols[i])
	polygon(x=c(nG_vals,rev(nG_vals)), 
		y=c(plotdat[[i]][1,],rev(plotdat[[i]][3,])),border=NA,col=fills[i])
}
dev.off()


