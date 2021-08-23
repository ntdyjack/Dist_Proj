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
results_rna <- lapply(gsets, function(x) sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,x,'.',y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  tmp <- tmp[1]/tmp[2]
  return(tmp) 
}))

results_sim <- lapply(gsets, function(x) sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,'splat.',x,'.',y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  tmp <- tmp[1]/tmp[2]
  return(tmp)
}))

#rna
xlim <- c(0.75,1.05) #c(0.0,3.0) 
xticks <- pretty(xlim,n=5)
xlabs <- formatC(xticks,format='f',digits=2)

#sim
ylim <- c(0.95,1.05)
yticks <- pretty(ylim,n=5) #c(0.00,1.00,2.00,3.00)
ylabs <- formatC(yticks,format='f',digits=2)


palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))
pch_labs <- rep(c(15,16,17),5)

pdf("/users/ndyjack/Dist_Proj/images/gapstat_rnaseq_vs_splatter.pdf",width=8,height=8)
for(j in 1:4){
  #dist legend (top)
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')


  #real vs splatter scatterplot
  par(new = "TRUE",plt = c(0.10,0.95,0.10,0.80))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim,xlim=xlim,axes=F,type='n',pch=16,col='black')
  abline(h=1,v=1,col='black',lty=2)
  points(x=results_rna[[j]],y=results_sim[[j]],pch=pch_labs,col=col_labs,cex=1.5)
  axis(side=2,labels=ylabs,at=yticks,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xlabs,at=xticks,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text="G(1)/G2 (simulation)",cex=1.2,line=2.5,font=2,las=3)
  mtext(side=1,text="G(1)/G2 (scRNA-seq)",cex=1.2,line=2.5,font=2,las=1)
  mtext(side=3,text=gsets[j],cex=1.2,line=0.2,font=2,las=1)
}
dev.off()

