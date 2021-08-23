#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)
rndup <- function(x,d) { ceiling(x*d) / d }
rnddn <- function(x,d) { floor(x*d) / d }

#trajcor
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
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5, 
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

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


results_dir <- '/users/ndyjack/Dist_Proj/tables/trajcor_new/'
#trajcor on real data
cor_sorted <- results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- read.table(paste0(results_dir,y,".",x,".txt"),header=T)
  return(mean(tmp$pearson))
}))
cor_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(results_sorted[grep(y,rownames(results_sorted)),x])
}))
colnames(cor_means) <- names(dists_vec)

res_dir <- "/users/ndyjack/Dist_Proj/tables/trajcor_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_tc <- sapply(res_files, function(x) {
  tmp <- read.table(paste0(res_dir,x),header=T)
  tmp <- mean(tmp$pearson)
  return(tmp)
})
pltdt_nmt <- unlist(lapply(1:3, function(i) lapply(1:3, function(j) paste0('t',i,'_','d',j,'_hvg1k.mtx.'))),recursive=T)
pltdt_tct <- lapply(pltdt_nmt, function(y) sapply(dists_vec, function(x) res_tc[grep(paste0(y,x,'\\.'),names(res_tc))]))

ylim_tc <- c(0.00,1.0)
yticks_tc <- c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_tc <- formatC(yticks_tc,format='f',digits=1)

palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))


#plotting parameters
zlim <- ylim_tc
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
ncol <- length(dists_vec)
nrow <- 4
xlim <- c(0,ncol+2)
ylim <-c(0,nrow+1)
#zleglocs <- c(0.01,0.25,0.50,0.75,0.99)
#zleglabs <- c("0.00","0.25","0.50","0.75","1.00")
#dist_gcol <- ifelse(dists_grp1=="Probabilistic",'red','blue')
#dist_pch <- ifelse(dists_grp1=="Probabilistic",15,16)

xlimbs <- c(0.3,6.7)
atl <- c(0.8,1.5,2.2, 2.8,3.5,4.2, 4.8,5.5,6.2)
atl[1:3] <- atl[1:3] - 0.3
atl[7:9] <- atl[7:9] + 0.3

pch_labs <- rep(c(15,16,17),5)

d_vec <- c('Easy','Medium','Hard')
t_vec <- c('Linear','Branched1','Branched2')

pdf("/users/ndyjack/Dist_Proj/images/trajcor_plot_v3.pdf",width=12,height=12)
  #dist legend (top)
  plot.new()
  par(new = "TRUE",plt = c(0.10,0.95,0.90,1.00),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  legend('top',legend=names(dists_vec),col=col_labs,ncol=5,
    x.intersp=0.5,y.intersp=1,text.width=0.15,cex=1.2,
    pch=pch_labs,pt.cex=2.0,bty='n')

  #beeswarm (middle)
  #trajcor x t
  par(new = "TRUE",plt = c(0.10,0.95,0.65,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlimbs,ylim=ylim_tc,xlab="",ylab="")
  beeswarm(pltdt_tct,at=atl,cex=1.5,pch=16,pwcol=rep(col_labs,9),pwpch=rep(pch_labs,9),corralWidth=0.4,corral='gutter',add=T)
  axis(side=2,labels=ylabs_tc,at=yticks_tc,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text=rep(d_vec,3),cex=0.8,line=0.1,font=2,at=atl,las=1,mgp=c(3, .5, 0))
  mtext(side=2,text='TrajCor',cex=1.2,line=2,font=2)
  mtext(side=1,text=t_vec,cex=1.2,line=1.2,font=2,at=atl[c(2,5,8)],las=1,mgp=c(3, .5, 0))

  #heatmap (bottom)
  tmp <- cor_means
  o <- order(colMeans(tmp),decreasing=F)
  tmp <- tmp[,o]
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.65),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=1-0.7,labels=colnames(tmp))
  text(x=0.2,y=(1:nrow)-0.25,labels=sub("_","",gsets),srt=90)
  #text(x=ncol+1.0,y=(1:nrow),labels=rownames(tmp),col=dist_gcol[o])  
dev.off()


