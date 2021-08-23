library(RColorBrewer)
library(beeswarm)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

n <- length(dists_vec)

#col_labs <- gg_color_hue(n)
palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, n)

pch_labs <- rep(c(15,16,17),5)

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

#simulation params
n_vec <- sapply(c(1000,5000,10000),function(x) paste0('n=',x))
k_vec <- sapply(c(2,5,10), function(x) paste0('k=',x))
d_vec <- c('Easy','Medium','Hard')
b_vec <- c('Uniform','Proportional','Unbalanced')
params_test <- c('n','k','b')

#read in simulation results
res_dir <- "/users/ndyjack/Dist_Proj/tables/ari_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_ari <- sapply(res_files, function(x) {
  tmp <- readLines(paste0(res_dir,x))
  tmp <- strsplit(tmp,' ')[[1]][2]
  tmp <- as.numeric(tmp)
  return(tmp)
})
pltdt_nmn <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('n',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_arin <- lapply(pltdt_nmn, function(y) sapply(dists_vec, function(x) res_ari[grep(paste0(y,x,'\\.'),names(res_ari))]))
pltdt_nmk <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('k',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_arik <- lapply(pltdt_nmk, function(y) sapply(dists_vec, function(x) res_ari[grep(paste0(y,x,'\\.'),names(res_ari))]))
pltdt_nmb <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('b',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_arib <- lapply(pltdt_nmb, function(y) sapply(dists_vec, function(x) res_ari[grep(paste0(y,x,'\\.'),names(res_ari))]))

ylim_ari <- c(-0.005,1.0)
yticks_ari <- c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_ari <- formatC(yticks_ari,format='f',digits=1)

#
res_dir <- "/users/ndyjack/Dist_Proj/tables/kacc_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_files <- grep('k50', res_files, value=T)
res_kacc <- sapply(res_files, function(x) {
  tmp <- readLines(paste0(res_dir,x))
  tmp <- strsplit(tmp,' ')[[1]][2]
  tmp <- as.numeric(tmp)
  return(tmp)
})
pltdt_nmn <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('n',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_kaccn <- lapply(pltdt_nmn, function(y) sapply(dists_vec, function(x) res_kacc[grep(paste0(y,x,'\\.'),names(res_kacc))]))
pltdt_nmk <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('k',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_kacck <- lapply(pltdt_nmk, function(y) sapply(dists_vec, function(x) res_kacc[grep(paste0(y,x,'\\.'),names(res_kacc))]))
pltdt_nmb <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('b',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_kaccb <- lapply(pltdt_nmb, function(y) sapply(dists_vec, function(x) res_kacc[grep(paste0(y,x,'\\.'),names(res_kacc))]))

ylim_kacc <- c(0.10,1.0)
yticks_kacc <- c(0.10,0.30,0.50,0.70,0.90)
ylabs_kacc <- formatC(yticks_kacc,format='f',digits=1)

#
res_dir <- "/users/ndyjack/Dist_Proj/tables/gplus_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_gplus <- sapply(res_files, function(x) {
  tmp <- readLines(paste0(res_dir,x))
  tmp <- strsplit(tmp,' ')[[1]][2]
  tmp <- 1 - as.numeric(tmp)
  return(tmp)
})
pltdt_nmn <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('n',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_gplusn <- lapply(pltdt_nmn, function(y) sapply(dists_vec, function(x) res_gplus[grep(paste0(y,x,'\\.'),names(res_gplus))]))
pltdt_nmk <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('k',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_gplusk <- lapply(pltdt_nmk, function(y) sapply(dists_vec, function(x) res_gplus[grep(paste0(y,x,'\\.'),names(res_gplus))]))
pltdt_nmb <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('b',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_gplusb <- lapply(pltdt_nmb, function(y) sapply(dists_vec, function(x) res_gplus[grep(paste0(y,x,'\\.'),names(res_gplus))]))

ylim_gplus <- c(0.70,1.00)
yticks_gplus <- seq(0.70,1.00,by=0.1)
ylabs_gplus <- formatC(yticks_gplus,format='f',digits=1)

xlim <- c(0.3,6.7)
atl <- c(0.8,1.5,2.2, 2.8,3.5,4.2, 4.8,5.5,6.2)
atl[1:3] <- atl[1:3] - 0.3
atl[7:9] <- atl[7:9] + 0.3

pdf("/users/ndyjack/Dist_Proj/images/splatsim_plots_3metrics.pdf",width=12,height=12)

  #plot 1, vary n
  #dist legend
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')

  #ari x n
  par(new = "TRUE",plt = c(0.10,0.95,0.65,0.90),las = 1,cex.axis = 1)
  plot.new() 
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_ari,xlab="",ylab="")
  beeswarm(pltdt_arin,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='ARI',cex=1.2,line=2,font=2)

  #kacc x n
  par(new = "TRUE",plt = c(0.10,0.95,0.35,0.60),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_kacc,xlab="",ylab="")
  beeswarm(pltdt_kaccn,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_kacc,at=yticks_kacc,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='kAcc',cex=1.2,line=2,font=2)

  #1-G+ x n
  par(new = "TRUE",plt = c(0.10,0.95,0.05,0.30),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_gplus,xlab="",ylab="")
  beeswarm(pltdt_gplusn,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_gplus,at=yticks_gplus,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=n_vec,cex=1.2,line=1.2,font=2,at=atl[c(2,5,8)],las=1,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='1-G+',cex=1.2,line=2,font=2)


  #plot 2, vary k
  #dist legend
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')

  #ari x k
  par(new = "TRUE",plt = c(0.10,0.95,0.65,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_ari,xlab="",ylab="")
  beeswarm(pltdt_arik,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='ARI',cex=1.2,line=2,font=2)

  #kacc x k
  par(new = "TRUE",plt = c(0.10,0.95,0.35,0.60),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_kacc,xlab="",ylab="")
  beeswarm(pltdt_kacck,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_kacc,at=yticks_kacc,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='kAcc',cex=1.2,line=2,font=2)

  #1-G+ x k
  par(new = "TRUE",plt = c(0.10,0.95,0.05,0.30),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_gplus,xlab="",ylab="")
  beeswarm(pltdt_gplusk,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_gplus,at=yticks_gplus,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=k_vec,cex=1.2,line=1.2,font=2,at=atl[c(2,5,8)],las=1,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='1-G+',cex=1.2,line=2,font=2)

 #plot 3, vary b
 #dist legend
 plot.new()
 par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
 plot.new()
 plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
 legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
   x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
   pch=pch_labs,pt.cex=2.0,bty='n')

 #ari x b
 par(new = "TRUE",plt = c(0.10,0.95,0.65,0.90),las = 1,cex.axis = 1)
 plot.new()
 plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_ari,xlab="",ylab="")
 beeswarm(pltdt_arib,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
 axis(side=2,labels=ylabs_ari,at=yticks_ari,cex=1.0,las=2,mgp=c(3, .5, 0))
 mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
 mtext(side=2,text='ARI',cex=1.2,line=2,font=2)

 #kacc x b
 par(new = "TRUE",plt = c(0.10,0.95,0.35,0.60),las = 1,cex.axis = 1)
 plot.new()
 plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_kacc,xlab="",ylab="")
 beeswarm(pltdt_kaccb,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
 axis(side=2,labels=ylabs_kacc,at=yticks_kacc,cex=1.0,las=2,mgp=c(3, .5, 0))
 mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
 mtext(side=2,text='kAcc',cex=1.2,line=2,font=2)

 #1-G+ x b
 par(new = "TRUE",plt = c(0.10,0.95,0.05,0.30),las = 1,cex.axis = 1)
 plot.new()
 plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim_gplus,xlab="",ylab="")
 beeswarm(pltdt_gplusb,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
 axis(side=2,labels=ylabs_gplus,at=yticks_gplus,cex=1.0,las=2,mgp=c(3, .5, 0))
 mtext(side=1,text=b_vec,cex=1.2,line=1.2,font=2,at=atl[c(2,5,8)],las=1,mgp=c(3, .5, 0))
 mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
 mtext(side=2,text='1-G+',cex=1.2,line=2,font=2)

dev.off()


