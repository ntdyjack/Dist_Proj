##
#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,heatmap3,beeswarm,doParallel)

dists_vec <- c(
  'L1','L2','SQRL2','JSM/PSM','JSD/PSD','MKL',
  'POISSON','HELLINGER','BHAT_METR','BHAT_DIST',
  'TVD','LLR','EMD','REV_MKL','REV_POISSON','UWLLR','OLLR'
)
idx_r <- 1:length(dists_vec)
idx_r <- idx_r[-which(dists_vec=="EMD")]
idx_c <- idx_r - 1

#location of data files
data_dir <- "/users/ndyjack/Dist_Proj/rdata/"
data_names <- list.files(data_dir)
data_names <- unique(sapply(data_names, function(x) strsplit(x,'\\.')[[1]][1]))
data_names <- data_names[1]
rdata_names <- lapply(data_names, function(x) paste(x,idx_c,"rda",sep="."))
#location of identity files
ident_dir <- "/users/ndyjack/Dist_Proj/tables/"
ident_vec <- sapply(data_names, function(x) paste0(x,"_labs.txt")) 
#location of output files
output_dir <- "/users/ndyjack/Dist_Proj/images/"

#calculate full g+
i <- 1
cluster_lab <- as.numeric(readLines(paste0(ident_dir,ident_vec[i])))
Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
Nz <- (Nt*(Nt-1))/2
#coerce to paired vector (IE, upper triangular cluster idenity matrix)
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
ind_vec <- apply(ind_vec,1,function(x) {
  tmp <- rep(0,length(x))
  tmp[which(x)] <- cluster_lab[which(x)]
  return(tmp)
})
ind_vec <- ind_vec[upper.tri(ind_vec)]
ind_intra <- which(ind_vec>0)
ind_inter <- which(ind_vec==0)
#read in results (don't want to calculate it agagin)
gplus_results <- read.table(paste0(ident_dir,"gplus_results.txt"),header=F,sep=" ")
gplus_results <- sapply(idx_c, function(x) gplus_results[grep(paste0("\\.",x,"\\."),gplus_results[,1]),2]) 

#approximate g+ with two methods
# order statistics, take representive points (from min to max)
# randomly sample that amount

set.seed(1234)
n_orders <- c(1e3,5e3,1e4,5e4)
Nz_tmp <- sapply(n_orders, function(x)  ((x+x)*(x+x-1))/2 )

sm_results <- lapply(1:16, function(j) {
  dist_vec <- readRDS(paste0(data_dir,rdata_names[[i]][j]))[[1]]
  dist_inter <- sort(dist_vec[ind_inter])
  dist_intra  <- sort(dist_vec[ind_intra])
  rm(dist_vec) #memory cleanup
  orders_inter <- lapply(n_orders, function(x) round(seq(1,length(ind_inter),length.out=x)))
  orders_intra <- lapply(n_orders, function(x) round(seq(1,length(ind_intra),length.out=x)))
  randos_inter <- lapply(n_orders, function(x) sample(1:length(ind_inter),size=x,replace=F))
  randos_intra <- lapply(n_orders, function(x) sample(1:length(ind_intra),size=x,replace=F))

  sm_tmp <- sapply(1:length(n_orders), function(k) {
    sm_ord <- sum(sapply(dist_intra[orders_intra[[k]]], function(x) sum(x > dist_inter[orders_inter[[k]]]))) / Nz_tmp[k]
    sm_ran <- sum(sapply(dist_intra[randos_intra[[k]]], function(x) sum(x > dist_inter[randos_inter[[k]]]))) / Nz_tmp[k]
    return(c(sm_ord,sm_ran))
  })
  return(sm_tmp)
})

xlab <- "Number of samples"
ylab <- "G+"
xlim=c(min(n_orders),max(n_orders))
xticklocs <- xticklabs <- n_orders
nticks <- 4

pdf(paste0(output_dir,"approx_gp_jurkat.pdf"),width=12,height=12)
par(mfrow=c(4,4))
for(i in 1:16){
  plot.new()
  par(new = "TRUE",plt = c(0.20,0.95,0.80,0.90),las = 1, cex.axis = 1)
  plot(x=0,y=0,xlim=c(0,1),ylim=c(0,1),type='n',bty='n',xaxt="n",yaxt="n",xlab="",ylab="",main="")
  legend('bottom',col=c("black","blue","red"),legend=c("Ordered","Random","Truth"),
    lty=c(2,2,1),pch=c(16,16,NA),horiz=T,bty='n',cex=1.0)
  mtext(side=3,text=dists_vec[idx_r[i]],cex=1.2,font=2,line=0.0,las=1)
    
  par(new = "TRUE",plt = c(0.20,0.95,0.20,0.80),las = 1, cex.axis = 1)
  ylim <- c(min(sm_results[[i]]),max(sm_results[[i]]))
  yticklocs <- pretty(ylim,nticks)
  yticklabs <- formatC(yticklocs,digits=3,format='f')
  plot(x=n_orders,y=sm_results[[i]][2,],xaxt="n",yaxt="n",xlab="",ylab="",main="",
    type='b',col="blue",pch=16,lty=2,lwd=1,xlim=xlim,ylim=ylim)
  points(x=n_orders,y=sm_results[[i]][1,],xaxt="n",yaxt="n",xlab="",ylab="",main="",
    type='b',col="black",pch=16,lty=2,lwd=1)
  abline(h=gplus_results[i],lty=1,col="red",lwd=1.1)
  axis(side=2,labels=yticklabs,at=yticklocs,cex=1.0,,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticklabs,at=xticklocs,cex=1.0,las=1,mgp=c(3, .5, 0),line=0)
  mtext(side=1,text=xlab,cex=1.1,line=1.8,font=2,las=1)
  mtext(side=2,text=ylab,cex=1.1,line=2.5,font=2,las=3)
}
dev.off()

