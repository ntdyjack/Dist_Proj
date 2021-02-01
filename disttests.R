#load accessory packages, set parameters
set.seed(17462)
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,heatmap3)

#location of data file
counts_loc <- paste0(args[1],args[3],"_expr.mtx") 
counts_mat <- readMM(counts_loc)
counts_mat <- t(as.matrix(counts_mat))
cluster_lab <- as.integer(readLines(paste0(args[1],args[3],"_labs.txt")))
#change the number of cells you want to test here
#commment out the lines to test the entire dataset
counts_mat <- counts_mat[1:1000,]
cluster_lab <- cluster_lab[1:1000]

#names of the distances being tested
dists_vec <- c(
  'L1','L2','SQRL2','JSM/PSM','JSD/PSD','MKL',
  'POISSON','HELLINGER','BHAT_METRIC','BHAT_DISTANCE',
  'TVD','LLR','EMD','REV_MKL','REV_POISSON','UWLLR','OLLR'
)
idx_test <- 1:length(dists_vec)

cid_vec <- unique(cluster_lab)
K <- length(cid_vec)
cid_list <- lapply(cid_vec, function(i) which(cluster_lab==i))
N <- nrow(counts_mat)
#Nt is the number of distinct pairs
Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
#Nw is the  total number of distinct pairs within each cluster
Nw <- sum(sapply(cid_list, function(i) length(i)*(length(i)-1)/2))
#Nb is the total numer of distinct pairs outside of each cluster
Nb <- Nt - Nw
nu_zer <- (Nt*(Nt-1))/2
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
ind_vec <- apply(ind_vec,1,function(x) {
  tmp <- rep(0,length(x))
  tmp[which(x)] <- cluster_lab[which(x)]
  return(tmp)
})
#vector for cluster identities between all unique pairs (0 means not in same clsuter)
ind_vec <- ind_vec[upper.tri(ind_vec)]

results_mat <- sapply(idx_test, function(i) {
  cat("running... ",dists_vec[i],"\n")

  #calculate and time distance
  start_time <- Sys.time()
  dist_vec <- Rfgc::dist_matrixdd(counts_mat, i)
  time_tmp <- difftime(time2=start_time,time1=Sys.time(),units="secs")
  #vector for distances between all unique pairs
  dist_vec <- dist_vec[upper.tri(dist_vec)]
  #s- number of times intra-clsuter distances are strictly greater than inter-cluster-distances
  sm <- sum(sapply(which(ind_vec>0), function(j) sum(dist_vec[j] > dist_vec[which(ind_vec==0)])))   
  
  return(c("Time (s)" = time_tmp,"Gplus" = sm / nu_zer))
})
colnames(results_mat) <- dists_vec[idx_test]

cat("plotting results...\n")
##color function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#plotting parameters
cols <- gg_color_hue(length(idx_test))
cols <- sample(cols,length(cols))
fills <- paste0(cols,"35")
shps <- rep(21:25,5)[idx_test]
shps <- sample(shps,length(shps))
ltys <- rep(1:5,5)[idx_test]
tick_min <- 0.05
res_names <- rownames(results_mat)
title_lab <- paste0(args[3],", n=", nrow(counts_mat))
##legend function
draw_legend <- function(){
  par(new = "TRUE",plt = c(0.15,0.95,0.03,0.23),las = 1, cex.axis = 1)
  plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n",
    yaxt="n",bty='n',xlim=c(0,1),ylim=c(0,1))
  legend('left',legend=dists_vec[idx_test],col=cols,ncol=4,
    x.intersp=0.6,y.intersp=1,text.width=0.25,cex=1.0,
    pch=shps,pt.bg="#00000050",pt.cex=1.0,bty='n')
}
#axis information for time
min_time <- round(min(results_mat[1,]),0)
max_time <- round(max(results_mat[1,]),0) 
ticks_time <- pretty(c(min_time,max_time),4)
labs_time <- formatC(ticks_time,format="f",digits=0)


pdf(paste0(args[2],args[3],"_disttest_scatterplot.pdf"),width=8,height=6)
for(i in 2:nrow(results_mat)){
  tmp <- results_mat[i,] 
  ymin <- (round(min(tmp),2)/tick_min)*tick_min
  ymax <- (round(max(tmp),2)/tick_min)*tick_min 
  yticks <- pretty(c(ymin,ymax),4) 
  ylabs <- formatC(yticks,format="f",digits=2)
  plot.new()
  par(new = "TRUE",plt = c(0.15,0.90,0.30,0.95))
  plot(x=results_mat[1,],y=tmp,
    ylim=c(min(yticks),max(yticks)),xlim=c(min(ticks_time),max(ticks_time)),
    main="",xlab="",ylab="",xaxt="n", yaxt="n",
    pch=shps,col=cols,bg="#00000050",cex=2)
  axis(side=1,labels=labs_time,at=ticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
  axis(side=2,labels=ylabs,at=yticks,cex=0.8,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=res_names[1],cex=1.0,line=1.5,font=2)
  mtext(side=2,text=res_names[i],cex=1.0,line=3.0,font=2,las=3)
  mtext(side=3,text=title_lab,cex=1.0,line=0.2,font=2)
  draw_legend()
}

heatmap3(results_mat[,order(results_mat[1,],decreasing=F)],
  scale="row",cexRow=1.2,cexCol=0.6,ColSideColors=cols,ColSideLabs="",Rowv=NA,Colv=NA)
dev.off()
