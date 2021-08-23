args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(RColorBrewer,viridis,leaflet)

#####
files = c(
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl",
  "cellbench10x_rnaseq_3cl",
  "cellbench10x_rnaseq_5cl",
  "cellbenchDS_rnaseq_3cl",
  "cellbenchCS2_rnaseq_3cl",
  "cellbenchCS2_rnaseq_5cl"
)
lf <- length(files)

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
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

results_dir <- "/users/ndyjack/Dist_Proj/tables/ari_new/"
res <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".kmeds.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
res_means <- sapply(1:length(dists_vec), function(x) sapply(sub("_","",gsets), function(y) {
  mean(res[grep(y,rownames(res)),x])
}))
ord <- order(colMeans(res_means),decreasing=F)
res <- res[,ord]

#plotting parameters
zlim <- c(min(res),max(res))
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
nc <- ncol(res)
nr <- nrow(res) 
xlim <- c(-0.5,nc+1)
ylim <-c(-0.5,nr+1)
zleglocs <- c(0.01,0.25,0.50,0.75,0.99)
zleglabs <- c("0.00","0.25","0.50","0.75","1.00")
txtcut <- 0.20
#y vertices for heatmap
yvals <- sapply(1:nr, function(i) c(i-0.5,i+0.5))
shrink <- 0.1
yvals2 <- lapply(1:lf, function(i) {
  tmp <- yvals[,files_clvec==i]
  tmp[1,] <- tmp[1,] - 0:(ncol(tmp)-1)*shrink
  tmp[2,] <- tmp[2,] - 1:ncol(tmp)*shrink
  return(tmp)
})
yvals <- do.call(cbind,yvals2)
ygenelocs <- unlist(lapply(yvals2, function(x) colMeans(x)))
ydatalocs <- unlist(lapply(yvals2, function(x) mean(x)))

pdf("/users/ndyjack/Dist_Proj/figures_final/supfig_5.1.pdf",width=8,height=12)

 #draw heatmap
  plot.new()
  par(new = "TRUE",plt = c(0.05,0.95,0.05,0.88),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:nc){
    for(j in 1:nr){
      rect(xleft=k-0.5,xright=k+0.5,ybot=yvals[1,j],ytop=yvals[2,j],border='black',col=col_pal_fxn(res[j,k]))
      text(x=k,y=mean(yvals[,j]),labels=formatC(res[j,k],digits=2,format='f'),cex=0.6,col=ifelse(res[j,k]>txtcut,'black','white'))
    }
  }
  text(x=1:nc,y=0.0,labels=colnames(res),cex=0.8)
  text(x=0.1,y=ygenelocs,labels=rep(sub("_","",gsets),lf),cex=0.5,srt=0)
  text(x=nc+0.7,y=ydatalocs,srt=90,cex=0.5,labels=gsub("_"," ",files))

  #color legend
  par(new = "TRUE",plt = c(0.20,0.80,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))
  mtext(side=3,text='ARI', cex=1.2)

dev.off()

write.csv(res,"/users/ndyjack/Dist_Proj/tables_final/suptab_7.1.csv",quote=F)
