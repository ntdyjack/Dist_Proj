map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}
library(plyr)

set.seed(1234)
n <- 1000
m <- 500
Nd <- (n*(n-1))/2
Nz <- (Nd*(Nd-1))/2
#((n+n)*(n+n-1))/2
pct_samps <- c(.000005,0.00005,.0005,.005,.05,.125) 
#seq(0.000025,0.001,by=0.000025)
#seq(0.00025,0.05,length.out=6) #c(.000025,.00025,.0025,.025)
#Nz_tmp <- sapply(n_orders, function(x)  ((x+x)*(x+x-1))/2 )
## 2 datasets (1 and 2 clusters)
#first dataset, assign random labels
dat1 <- sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1))
labs1 <- round(runif(n=n))
ind1 <- sapply(labs1, function(x) sapply(labs1, function(y) x==y))
#tmp <- ind1[upper.tri(ind1)]
ind1 <- ind1[upper.tri(ind1)]
ind1_intra <- which(ind1)
ind1_inter <- which(!ind1)
dist1 <- as.matrix(dist(t(dat1)))
dist1 <- dist1[upper.tri(dist1)]
intra1 <- dist1[ind1_intra]
inter1 <- dist1[ind1_inter]
#calculate true g+
ptm <- proc.time()
gp1_true <- 1-sum(sapply(intra1, function(x) sum(x>inter1))) / Nz
t1 <- proc.time() - ptm
t1 <- t1[3]

#estimate g+ using order statisttcs
#sort for estimator
intra1 <- sort(intra1)
inter1 <- sort(inter1)
#loop over perctange of distances to sample
gp1_estm <- lapply(pct_samps, function(x) {
  n_samp <- round(x*length(ind1))
  Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
  orders_inter <- round(seq(1,length(inter1),length.out=n_samp))
  orders_intra <- round(seq(1,length(intra1),length.out=n_samp))
  ptm <- proc.time()
  gp_samp <- 1-sum(sapply(intra1[orders_intra], function(x) sum(x > inter1[orders_inter]))) / Nz_samp
  t1 <- proc.time() - ptm  
  return(list(gp_samp,t1[3]))
})

#second dataset, real labels (based on simulated cluster)
dat2 <- cbind(
  sapply(1:(n/2), function(i) rnorm(n=m,mean=0.3,sd=1)),
  sapply(1:(n/2), function(i) rnorm(n=m,mean=-0.3,sd=1))
)
labs2 <- c(rep(0,n/2),rep(1,n/2))
ind2 <- sapply(labs2, function(x) sapply(labs2, function(y) x==y))
#tmp <- ind2[upper.tri(ind2)]
ind2 <- ind2[upper.tri(ind2)]
ind2_intra <- which(ind2)
ind2_inter <- which(!ind2)
dist2 <- as.matrix(dist(t(dat2)))
dist2 <- dist2[upper.tri(dist2)]
intra2 <- dist2[ind2_intra]
inter2 <- dist2[ind2_inter]
#calculate true g+
ptm <- proc.time()
gp2_true <- 1-sum(sapply(intra2, function(x) sum(x>inter2))) / Nz
t2 <- proc.time() - ptm
t2 <- t2[3]

#estimate g+ using order statisttcs
#sort for estimator
intra2 <- sort(intra2)
inter2 <- sort(inter2)
#loop over perctange of distances to sample
gp2_estm <- lapply(pct_samps, function(x) {
  n_samp <- round(x*length(ind2))
  Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
  orders_inter <- round(seq(1,length(inter2),length.out=n_samp))
  orders_intra <- round(seq(1,length(intra2),length.out=n_samp))
  ptm <- proc.time()
  gp_samp <- 1-sum(sapply(intra2[orders_intra], function(x) sum(x > inter2[orders_inter]))) / Nz_samp
  t2 <- proc.time() - ptm
  return(list(gp_samp,t2[3]))
})

pca1 <- prcomp(t(dat1))$x
pca2 <- prcomp(t(dat2))$x

#final plot, 4 rows x 2 columns
#row 1, PCA plots of 2 simulated datasets
#row 2, histgrams of inter and intra-cluster distances
#row 3, asymptotic G+ of estimator
#row 4, asymptotic time properties of estimator

cl1 <- '#0000ff64' #blue
cl2 <- '#ff000064' #red

cols1 <- ifelse(labs1==1,cl1,cl2)
cols2 <- ifelse(labs2==1,cl1,cl2)

bins <- seq(floor(min(c(dist1,dist2))),ceiling(max(c(dist1,dist2))),length.out=20)

plotlocs <- rbind(
  c(0.10,0.50,0.75,0.98), #row1 col1
  c(0.55,0.95,0.75,0.98), #row1 col2
  c(0.10,0.50,0.55,0.70), #row2 col1
  c(0.55,0.95,0.55,0.70), #row2 col2
  c(0.10,0.50,0.30,0.50), #row3 col1
  c(0.55,0.95,0.30,0.50), #row3 col2
  c(0.10,0.50,0.05,0.25), #row4 col1
  c(0.55,0.95,0.05,0.25)  #row4 col2
)

#gp1_true <- gp2_true <- 0.75

xlab_row1 <- 'PC1'
ylab_row1 <- 'PC2'
xlab_row2 <- 'Dist. (L2)'
ylab_row2 <- '% Dists'
xlab_row3 <- '% Dists. Sampled'
ylab_row3 <- '1-G+'
ylab_row4 <- 'Time (s)'

yticks_row2 <- c(0.0,0.10,0.20,0.30)
ytclbs_row2 <- formatC(yticks_row2,digits=2,format='f')

xticks_row2 <- pretty(bins,5)

xlims_row3 <- c(-.05,1.05)
xticks_row3 <- c(0.0,0.15,0.30,0.50,0.75,1.0) #pct_samps*200
xtclbs_row3 <- c('.0001', '0.01', '0.1','1','10','25')



ylims_row3 <- c(0.65,1.00)
yticks_row3 <- c(0.70,0.80,0.80,0.90,1.00)
ytclbs_row3 <- formatC(yticks_row3,digits=1,format='f')


et1 <- sapply(gp1_estm, function(x) x[[2]])
et2 <- sapply(gp2_estm, function(x) x[[2]])

tmin <- round_any(min(c(et1[5:6],et2[5:6])),10,f=floor)
#round_any(min(c(t1,t2)),10,f=floor)
tmax <- round_any(max(c(t1,t2)),10,f=ceiling)

et1[5:6] <- map_fxn(et1[5:6],0.8,2,tmin,tmax)
et2[5:6] <- map_fxn(et2[5:6],0.8,2,tmin,tmax)

t1_scld <- map_fxn(t1,0.8,2,tmin,tmax)
t2_scld <- map_fxn(t2,0.8,2,tmin,tmax)
prtbig <- round_any(seq(tmin,tmax,length.out=3),10)
#seq(tmin,tmax,length.out=3) #pretty(tmin,tmax,4)

ylims_row4 <- c(-0.1,2.1) #c(tmin-2.1,tmax+2.1)
yticks_row4 <- c(0.0,0.5,seq(from=0.8,to=2.0,length.out=length(prtbig)))
ytclbs_row4 <- c('0.0','0.5',as.character(prtbig))
#c(tmin-2,tmin-1.5,tmin-1.0,pretty(ylims_row4,2))
#ytclbs_row4 <- c('0.0','0.5','1.0',as.character(yticks_row4[-c(1:3)])) 
#et1 <- sapply(gp1_estm, function(x) x[[2]]) + 
#et2 <- sapply(gp2_estm, function(x) x[[2]]) + 
 #10^ceiling(log10(mean(c(t1,t2)*.05)))
#lt2s <- log10(sapply(gp2_estm, function(x) x[[2]])+1.01)
#t1 <- t2 <- 680
#ylims_row4 <- #c(600,700)
#yticks_row4 <- c(600,600.5,601,650,700)


#ylims_row4 <- c(-0.1,5)
#yticks_row4s <- c(0,50,300,700)
#yticks_row4 <- log10(yticks_row4s+1.01)  #c(0,200,400,600) #pretty(ylims_row4,4)
#ytclbs_row4 <- as.character(yticks_row4s) 

pdf('/users/ndyjack/Dist_Proj/figures_final/gplus_supfig.pdf',height=10,width=7)
  plot.new()
  par(new = "TRUE",plt = plotlocs[1,],las = 1,cex.axis = 1)
  plot(x=pca1[,1],y=pca1[,2],pch=16,col=cols1,cex=0.7,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  mtext(side=1,text=xlab_row1,line=0.5,las=1,font=2)
  mtext(side=2,text=ylab_row1,line=0.5,las=3,font=2)

  par(new = "TRUE",plt = plotlocs[2,],las = 1,cex.axis = 1) 
  plot(x=pca2[,1],y=pca2[,2],pch=16,col=cols2,cex=0.7,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  mtext(side=1,text=xlab_row1,line=0.5,las=1,font=2)

  par(new = "TRUE",plt = plotlocs[3,],las = 1,cex.axis = 1)
  hist(x=intra1,breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',col=cl1,
    freq=F,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n')
  hist(x=inter1,breaks=bins,add=T,border='red',col=cl2,freq=F)
  mtext(side=1,text=xlab_row2,line=1.2,las=1,font=2)
  mtext(side=2,text=ylab_row2,line=2.5,las=3,font=2) 
  axis(side=1,at=xticks_row2,las=1,mgp=c(3, .3, 0),line=0.1)
  axis(side=2,labels=ytclbs_row2,at=yticks_row2,cex=1.0,las=2,mgp=c(3, .5, 0))

  par(new = "TRUE",plt = plotlocs[4,],las = 1,cex.axis = 1)
  hist(x=intra2,breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',col=cl1,
    freq=F,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n')
  hist(x=inter2,breaks=bins,add=T,border='red',col=cl2,freq=F)
  mtext(side=1,text=xlab_row2,line=1.2,las=1,font=2)
  axis(side=1,at=xticks_row2,las=1,mgp=c(3, .3, 0),line=0.1)

  par(new = "TRUE",plt = plotlocs[5,],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims_row3,ylim=ylims_row3,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks_row3, y=sapply(gp1_estm, function(x) x[[1]]),pch=1,col='blue')
  abline(h=gp1_true,lty=2,col="black",lwd=1.1)
  axis(side=2,labels=ytclbs_row3,at=yticks_row3,cex=1.0,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xtclbs_row3[1:3],at=xticks_row3[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[4:6],at=xticks_row3[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  #axis(side=1,labels=xtclbs_row3,at=xticks_row3,cex.axis=0.5,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab_row3,line=1.0,font=2,las=1)
  mtext(side=2,text=ylab_row3,line=2.5,font=2,las=3)

  par(new = "TRUE",plt = plotlocs[6,],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims_row3,ylim=ylims_row3,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks_row3, y=sapply(gp2_estm, function(x) x[[1]]),pch=1,col='blue')
  abline(h=gp2_true,lty=2,col="black",lwd=1.1)
  #axis(side=1,labels=xtclbs_row3,at=xticks_row3,cex=1.0,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[1:3],at=xticks_row3[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[4:6],at=xticks_row3[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab_row3,line=1.0,font=2,las=1)

  par(new = "TRUE",plt = plotlocs[7,],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims_row3,ylim=ylims_row4,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks_row3, y=et1,pch=1,col='blue')
  abline(h=t1_scld,lty=2,col="black",lwd=1.1)
  axis(side=2,labels=ytclbs_row2[1:2],at=yticks_row4[1:2],cex=1.0,las=2,mgp=c(3, .5, 0))
  axis(side=2,labels=ytclbs_row4[3:5],at=yticks_row4[3:5],cex.axis=1.0,las=2,mgp=c(3, .5, 0))
  #axis(side=1,labels=xtclbs_row3,at=xticks_row3,cex.axis=0.5,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[1:3],at=xticks_row3[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[4:6],at=xticks_row3[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab_row3,line=1.0,font=2,las=1)
  mtext(side=2,text=ylab_row4,line=2.5,font=2,las=3)

  par(new = "TRUE",plt = plotlocs[8,],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims_row3,ylim=ylims_row4,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks_row3, y=et2,pch=1,col='blue')
  abline(h=t2_scld,lty=2,col="black",lwd=1.1)
  #axis(side=1,labels=xtclbs_row3,at=xticks_row3,cex.axis=0.5,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[1:3],at=xticks_row3[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs_row3[4:6],at=xticks_row3[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab_row3,line=1.0,font=2,las=1)
dev.off()


