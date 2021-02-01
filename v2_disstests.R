##color function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
##
##legend function
draw_legend <- function(points=T){
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.23),las = 1, cex.axis = 1)
  plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n", yaxt="n",bty='n')
  if(points){
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=4,
      x.intersp=0.6,y.intersp=1,text.width=0.2,cex=1.0,
      pch=shps,pt.bg="#00000050",pt.cex=1.0,bty='n')
  }else{
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=4,
      x.intersp=0.5,y.intersp=1,text.width=0.2,cex=1.0,
      lty=ltys,lwd=1.5,bty='n')
  }
}

#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(Matrix,Rfgc,heatmap3)

#location of data file
#cluster_lab <-  as.numeric(readLines("/users/ndyjack/Dist_Proj/tables/jurkat293t_2cl_labs.txt"))
#counts_loc <- "/users/ndyjack/Dist_Proj/tables/jurkat293t_2cl_expr.mtx"
counts_loc <- paste0(args[1],args[3],"_expr.mtx") 
counts_mat <- readMM(counts_loc)
counts_mat <- t(as.matrix(counts_mat))
cluster_lab <- as.integer(readLines(paste0(args[1],args[3],"_labs.txt")))
counts_mat <- counts_mat[1:100,]
cluster_lab <- cluster_lab[1:100]

#names of the distances being tested
dists_vec <- c('L1','L2','SQRL2','JSM/PSM','JSD/PSD','MKL','POISSON','HELLINGER',
  'BHAT_METRIC','BHAT_DISTANCE','TVD','LLR',
  'EMD','REV_MKL','REV_POISSON','UWLLR','OLLR')
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
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
ind_vec <- apply(ind_vec,1,function(x) {
  tmp <- rep(0,length(x))
  tmp[which(x)] <- cluster_lab[which(x)]
  return(tmp)
}) 
ind_vec <- ind_vec[upper.tri(ind_vec)]
nu_zer <- (Nt*(Nt-1))/2
#nu_two <- ((Nb*(Nb-1))/2) + ((Nw*(Nw-1))/2) 

results_list <- lapply(idx_test, function(i) {
  cat("running... ",dists_vec[i],"\n")

  #calculate and time distance
  start_time <- Sys.time()
  dist_vec <- Rfgc::dist_matrixdd(counts_mat, i)
  time_tmp <- Sys.time() - start_time
  dist_vec <- dist_vec[upper.tri(dist_vec)]

  #s+, number of times intra-clsuter distances are strictly smaller than inter-cluster-distances
  sp <- sum(sapply(which(ind_vec>0), function(j) sum(dist_vec[j] < dist_vec[which(ind_vec==0)])))
  #s- number of times intra-clsuter distances are strictly greater than inter-cluster-distances
  sm <- sum(sapply(which(ind_vec>0), function(j) sum(dist_vec[j] > dist_vec[which(ind_vec==0)])))   
  #sum of number of intra-inter cluster ties within each group (times itself / 2)
  #nu_one <- sum(sapply(1:K, function(j) {
  #  tmp_idx <- which(ind_vec==j)
  #  ti <- sum(sapply(tmp_idx, function(l) sum(dist_vec[l] == dist_vec[which(ind_vec==0)])))
  #  return( ((ti*(ti-1))/2) )
  #}))

  return(c(
    "Time (s)" = time_tmp,
    "Gamma" = (sp - sm) / (sp + sm),
    "Gplus" = sm / nu_zer#,
#    "Tau" = (sp - sm) / sqrt((nu_zer - nu_one) * (nu_zer - nu_two))
  ))
})

cat("plotting results...\n")
#plotting parameters
cols <- gg_color_hue(length(idx_test))
fills <- paste0(cols,"35")
shps <- rep(21:25,5)[idx_test]
ltys <- rep(1:5,5)[idx_test]
tick_min <- 0.05
index_length <- length(results_list[[1]])
index_names <- names(results_list[[1]])
title_lab <- paste0(args[3],", n=", nrow(counts_mat))



pdf(paste0(args[2],args[3],"_disttest_v4.pdf"),width=8,height=6)
for(i in 1:index_length){
  tmp <- sapply(results_list, function(x) x[i])
#  tmp <- tmp/max(tmp)
  ymin <- (round(min(tmp),2)/tick_min)*tick_min
  ymax <- (round(max(tmp),2)/tick_min)*tick_min 
  yticks <- pretty(c(ymin,ymax),4) 
  ylabs <- formatC(yticks,format="f",digits=2)
  
  plot.new()
  par(new = "TRUE",plt = c(0.13,0.95,0.30,0.98))
  plot(x=1:length(tmp),y=tmp,ylim=c(min(yticks),max(yticks)),
    main="",xlab="",ylab="",xaxt="n", yaxt="n",
    pch=shps,col=cols,bg="#00000050",cex=2)
#  axis(side=1,labels=ylabs_time,at=yticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
  axis(side=2,labels=ylabs,at=yticks,cex=0.8,las=2,mgp=c(3, .5, 0))
#  mtext(side=1,text=yaxis_time,cex=1.0,line=1.5,font=2)
  mtext(side=2,text=index_names[i],cex=1.0,line=3.5,font=2,las=3)
  mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
  draw_legend(T)
}

tmp <- do.call(rbind,results_list)
rownames(tmp) <- dists_vec[idx_test]
heatmap3(t(tmp),scale="row",cexRow=0.6,cexCol=0.6,ColSideColors=cols)
dev.off()
