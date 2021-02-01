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

#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,randomcoloR)

#####
#read dist names from Rfcg output
#read input directory
#read in data files
#compute metric
#plot (beeswarm)
#datasets

files = c(
  "cellmix_rnaseq_traj",
  "rnamixcs2_rnaseq_traj",
  "rnamixss_rnaseq_traj"
)

gsets = c(
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
)
gsets <- rev(gsets)

datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQRL2'=2, 'JSM/PSM'=3,
  'JSD/PSD'=4, 'MKL'=5, 'POISSON'=6, 'HELLINGER'=7,
  'BHAT_METR'=8, 'BHAT_DIST'=9, 'TVD'=10,'LLR'=11,
  'REV_MKL'=12, 'REV_POIS'=13 , 'ITA_SAI'=16,
  'REV_ITA_SAI'=17, 'COS_DIST'=18, 
  #'COS_SIM'=20,'PRB_COS_SIM'=21,
  'PL2'=28,'PSL2'=29,
  'std_pca'='std_pca','trn_pca'='trn_pca','glm_pca'='glm_pca'
)

#location of results files
results_dir <- "/users/ndyjack/Dist_Proj/tables/trajcor/"
#result_tab <- read.table(result_, sep=" ",header=F)
#result_tab[,1] <- sub("\\.full","",result_tab[,1])
results_sorted <- lapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- read.table(paste0(results_dir,y,".",x,".txt"),header=T)
  return(tmp$pearson)
}))

ylim <- c(0.00,1.00)
yticklocs <-  pretty(ylim,4)
yticklabs  <-  as.character(formatC(yticklocs,digits=1,format='f'))
set.seed(1)
cols <- distinctColorPalette(k=length(dists_vec))
fills <- paste0(cols,"35")
xlim <- c(0,round(length(dists_vec)*1.5))
atl <- seq(from=xlim[1]+0.5,to=xlim[2]-0.5,length.out=length(dists_vec))

pdf("/users/ndyjack/Dist_Proj/images/trajcor_plots.pdf",width=10,height=6,family='sans')
for(x in gsets){
  tmp <- lapply(results_sorted, function(y) y[,grep(x,colnames(y))])
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.99,0.20,0.95))
  beeswarm(tmp,xaxt="n",yaxt="n",xlab="",ylab="",main="",pch=16,bg=fills,col=cols,ylim=ylim,at=atl,xlim=xlim)
  bxplot(tmp,add=T,at=atl)
  axis(side=2,labels=yticklabs,at=yticklocs,cex=1.5,las=2,mgp=c(3, .5, 0),font=1)
  mtext(side=2,text="Pearson Correlation",cex=1.2,line=2.5,font=2,las=3)
  mtext(side=1,text=names(dists_vec),cex=0.9,line=0.5,font=1,las=3,at=atl)
  mtext(side=3,text=sub("_","",x),cex=2.0,line=0.1)
}
dev.off()

