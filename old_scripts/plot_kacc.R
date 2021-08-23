##legend function
draw_legend <- function(points=T){
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.23),las = 1, cex.axis = 1)
  plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n", yaxt="n",bty='n')
  if(points){
    legend('left',legend=names(dists_vec),col=cols,ncol=4,
      x.intersp=0.6,y.intersp=1,text.width=0.17,cex=1.0,
      pch=shps,pt.bg="#00000050",pt.cex=1.0,bty='n',lty=ltys)
  }else{
    legend('left',legend=names(dists_vec),col=cols,ncol=4,
      x.intersp=0.5,y.intersp=1,text.width=0.2,cex=1.0,
      lty=ltys,lwd=1.5,bty='n')
  }
}

set.seed(2185)

#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer)

#####
#read dist names from Rfcg output
#read input directory
#read in data files
#compute metric
#plot (beeswarm)
#datasets


files = c(
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl",
  "cellbench10x_rnaseq_3cl",
  "cellbench10x_rnaseq_5cl",
  "cellbenchDS_rnaseq_3cl",
  "cellbenchCS2_rnaseq_3cl",
  "cellbenchCS2_rnaseq_5cl"
)

gsets = c(
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
)
gsets <- rev(gsets)

files_clvec <- as.vector(sapply(1:length(files), function(i) rep(i,4)))

datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQRL2'=2, 'JSM/PSM'=3,
  'JSD/PSD'=4, 'MKL'=5, 'POISSON'=6, 'HELLINGER'=7,
  'BHAT_METR'=8, 'BHAT_DIST'=9, 'TVD'=10,'LLR'=11,
  #'REV_MKL'=12, 'REV_POIS'=13 , 
  'ITA_SAI'=16,
  #'REV_ITA_SAI'=17, 
  'COS_DIST'=18, 
  #'COS_SIM'=20,'PRB_COS_SIM'=21,
  'PL2'=28,'PSL2'=29#,
#  'std_pca'='std_pca','trn_pca'='trn_pca','glm_pca'='glm_pca'
)

#location of results files
results_dir <- "/users/ndyjack/Dist_Proj/tables/kacc/"
#result_tab <- read.table(result_, sep=" ",header=F)
#result_tab[,1] <- sub("\\.full","",result_tab[,1])
results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".k30.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
idx <- seq(2,nrow(results_sorted),by=4)
heatmap_labs <- rep("",nrow(results_sorted))
heatmap_labs[idx] <-  sapply(rownames(results_sorted)[idx], function(x) paste0(strsplit(x,"_")[[1]][1:3],collapse="\n"))
#rownames(results_sorted)[idx] <- sapply(rownames(results_sorted)[idx], function(x) paste0(paste0(strsplit(x,"_")[[1]][1:3],collapse="_"),"\n",strsplit(x,"_")[[1]][4]))
#rownames(results_sorted)[-idx] <- sapply(rownames(results_sorted)[-idx], function(x) strsplit(x,"_")[[1]][4])
#rownames(results_sorted) <- gsub("rnaseq_|citeseq_|[0-9]cl_","",rownames(results_sorted))
#rownames(results_sorted)[!grepl("hvg",rownames(results_sorted))] <- paste0(rownames(results_sorted)[!grepl("hvg",rownames(results_sorted))],"_full")
#rownames(results_sorted) <- sub("_","\n",rownames(results_sorted))
#rownames(results_sorted) <- gsub("_","",rownames(results_sorted))

dist_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))
ylim <- c(0.60,1.00)
yticklocs <-  pretty(ylim,4)
yticklabs  <-  as.character(yticklocs)
xlim <- c(1,4)
#xlabs <- rownames(dist_means) #c("",
xticklocs <- 1:4
xticklabs <- rownames(dist_means)
palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
cols <- sample(palette3_all, length(dists_vec))
fills <- paste0(cols,"35")
shps <- rep(21:25,5)[1:length(dists_vec)]
ltys <- rep(1:5,5)[1:length(dists_vec)]

pdf("/users/ndyjack/Dist_Proj/images/kacc_plots.pdf")
heatmap3(t(results_sorted),scale="none",cexRow=0.8,cexCol=1,lasCol=1,Colv=NA,
  labCol=heatmap_labs,RowSideColors=cols,RowSideLabs='',ColSideColors=files_clvec,ColSideLabs='')
#mtext(at=seq(0,1,length.out=50),seq(0,1,length.out=50))
#spaghetti plot

  plot.new()
  par(new = "TRUE",plt = c(0.13,0.95,0.30,0.98))
  plot(x=1,y=1,xaxt="n",yaxt="n",xlab="",ylab="",main="",
     ylim=ylim,xlim=xlim,axes=F,type='n',pch=16,col='black',cex=1.1)
  for(i in 1:length(dists_vec)){
    points(x=1:4,y=dist_means[,i],type='b',pch=shps[i],lty=ltys[i],col=cols[i],bg=fills[i])
  }
  axis(side=2,labels=yticklabs,at=yticklocs,cex=1.1,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xticklabs,at=xticklocs,cex=1.1,las=1,mgp=c(3, .5, 0))
  mtext(side=1,text="# Genes",cex=1.2,line=1.5,font=2,las=1)
  mtext(side=2,text="Mean kacc (9 datasets)",cex=1.2,line=2.5,font=2,las=3)

  draw_legend(T)
dev.off()













#cat("plotting results...\n")
#plotting parameters
#cols <- gg_color_hue(length(dists_vec))
#fills <- paste0(cols,"35")
#shps <- rep(21:25,5)
#ltys <- rep(1:5,5)

#tick_min <- 0.05
#ymin <- (round(min(output_metrics),2)/tick_min)*tick_min
#ymax <- (round(max(output_metrics),2)/tick_min)*tick_min
#yticks <- pretty(c(ymin,ymax),4)
#ylabs <- formatC(yticks,format="f",digits=2)

#pdf(paste0(output_dir,"beeswarm_disttest.pdf"),width=8,height=6)
#beeswarm(output_metrics,pch=ggdat_pch,at=atl,add=T,cex=1.0,corral="omit",corralWidth=0.5)
#  axis(side=2,labels=formatC(ytxl,digits=1,format="f"),cex=0.7,line=0,at=ytxl,las=2,mgp=c(0.0,0.6,0))
#  mtext(text="",side=2,cex=1.5,line=1.8)
#  mtext(text=final_ref,side=1,cex=1.5,line=0.1,at=atl,mgp=c(0.0,0.4,0))
  
#  plot.new()
#  par(new = "TRUE",plt = c(0.13,0.95,0.30,0.98))
#  plot(x=1:length(tmp),y=tmp,ylim=c(min(yticks),max(yticks)),
#    main="",xlab="",ylab="",xaxt="n", yaxt="n",
#    pch=shps,col=cols,bg="#00000050",cex=2)
#  axis(side=1,labels=ylabs_time,at=yticks_time,cex=0.8,las=1,mgp=c(3, .5, 0))
#  axis(side=2,labels=ylabs,at=yticks,cex=0.8,las=2,mgp=c(3, .5, 0))
#  mtext(side=1,text=yaxis_time,cex=1.0,line=1.5,font=2)
#  mtext(side=2,text=index_names[i],cex=1.0,line=3.5,font=2,las=3)
#  mtext(side=3,text=title_lab,cex=1.0,line=0.3,font=2)
#  draw_legend(T)
#}

#tmp <- do.call(rbind,results_list)
#rownames(tmp) <- dists_vec[idx_test]
#heatmap3(t(tmp),scale="row",cexRow=0.6,cexCol=0.6,ColSideColors=cols)
#dev.off()

