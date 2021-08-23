#read in simulation results
library(Rtsne)
in_dir <- "/users/ndyjack/Dist_Proj/rdata/"
sim_files <- list.files(in_dir)
sim_files <- sim_files[grep("splatterTests.dr",sim_files)]
sim_files <- paste0(in_dir, sim_files)
sim_results <- readRDS(sim_files)
sim_results <- unlist(sim_results,recursive=F)

#plotting params
titles <- c("Euclidean Distance","Multinomial Distance","Poisson Distance")
seed <- 1427
nticks <- 4

set.seed(seed)
pdf("/users/ndyjack/Dist_Proj/images/splatterSim_tSNE.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4,4,4))
for(i in 1:length(sim_results)){
  col <- ifelse(as.character(sim_results[[i]]$id)=="Group1","blue","red")
  info <- unlist(strsplit(sim_results[[i]][[5]], ", "))
  info <- c(paste0(info[1:2],collapse=", "),
    paste0(info[3:4],collapse=", "),paste0(info[5:7],collapse=", "))
  for(j in 1:3){
    #plot results / params
    tsne_results <- Rtsne(sim_results[[i]][[j]],is_distance=T,pca=F)$Y
    xmin <- min(tsne_results[,1])
    xmax <- max(tsne_results[,1])
    ymin <- min(tsne_results[,2])
    ymax <- max(tsne_results[,2])
    xtickloc <- seq(xmin,xmax,length.out=nticks)
    ytickloc <- seq(ymin,ymax,length.out=nticks)
    xticklab <- formatC(xtickloc,format="f",digits=1)
    yticklab <- formatC(ytickloc,format="f",digits=1)

    plot(x=tsne_results[,1],y=tsne_results[,2],pch=16,cex=0.5,
      ylab="",yaxt="n",xlab="",xaxt="n",main="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),col=col)
    mtext(side=3,text=titles[j],font=2,line=0.1,cex=1.5)
    axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
    axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
    mtext(side=1,text=info[j],cex=0.8,line=2)
  }
}
dev.off()


pdf("/users/ndyjack/Dist_Proj/images/splatterSim_MDS.pdf",width=12,height=4)
par(mfrow=c(1,3),mar=c(4,4,4,4))
for(i in 1:length(sim_results)){
  col <- ifelse(as.character(sim_results[[i]]$id)=="Group1","blue","red")
  info <- unlist(strsplit(sim_results[[i]][[5]], ", "))
  info <- c(paste0(info[1:2],collapse=", "),
    paste0(info[3:4],collapse=", "),paste0(info[5:7],collapse=", "))
  for(j in 1:3){
    #plot results / params
    mds_results <- cmdscale(sim_results[[i]][[j]])
    xmin <- min(mds_results[,1])
    xmax <- max(mds_results[,1])
    ymin <- min(mds_results[,2])
    ymax <- max(mds_results[,2])
    xtickloc <- seq(xmin,xmax,length.out=nticks)
    ytickloc <- seq(ymin,ymax,length.out=nticks)
    xticklab <- formatC(xtickloc,format="f",digits=1)
    yticklab <- formatC(ytickloc,format="f",digits=1)

    plot(x=mds_results[,1],y=mds_results[,2],pch=16,cex=0.5,
      ylab="",yaxt="n",xlab="",xaxt="n",main="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),col=col)
    mtext(side=3,text=titles[j],font=2,line=0.1,cex=1.5)
    axis(side=1,labels=xticklab,at=xtickloc,cex=0.7,las=1,mgp=c(3, .5, 0))
    axis(side=2,labels=yticklab,at=ytickloc,cex=0.7,las=2,mgp=c(3, .7, 0))
    mtext(side=1,text=info[j],cex=0.8,line=2)
  }
}
dev.off()

#set up umap params (not a thing just yet)
seed <- 1427
library(umap)
umap_params <- umap.defaults
umap_params$metric <- "precomputed"
umap_params$random_state <- seed
umap_test <- umap(as.matrix(sim_results[[1]]$euc),umap_params)
