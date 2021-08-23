args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(RColorBrewer,viridis,leaflet)

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5,
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)

#simulation params
n_vec <- c(1000,5000,10000)
k_vec <- c(2,5,10)
d_vec <- c('Easy','Medium','Hard')
b_vec <- c('Unif.','Prop.','Unbal.')
params_test <- c('n','k','b')
datasets_vec <- unlist(lapply(params_test, function(x) sapply(1:3, function(i) sapply(1:3, function(j) paste0(x,j,'_d',i)))))
dataset_nmz <- c(
  rep(sapply(1:3, function(i) paste0('n=',n_vec[i])),3),
  rep(sapply(1:3, function(i) paste0('k=',k_vec[i])),3),
  rep(sapply(1:3, function(i) paste0('b=',b_vec[i])),3)
)
lf <- 9
files_clvec <- as.vector(sapply(1:lf, function(i) rep(i,3)))
#dataset_nmz <- as.vector(sapply(dataset_nmz, function(x) rep(x,3)))
#sapply(params_test, function(x) sapply(1:3, function(i) paste0(x,'=',n_vec[i])

results_dir <- "/users/ndyjack/Dist_Proj/tables/gplus_new/"
res <- 1 - sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,'splatsim_',y,'_hvg1k.mtx.',x,".txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
ord <- order(colMeans(res),decreasing=F)
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
ylim <-c(-0.5,nr+0.3)
zleglocs <- c(0.01,0.25,0.50,0.75,0.99)
zleglabs <- c("0.75","0.80","0.85","0.95","1.00")
txtcut <- 0.82
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
yparlocs <- unlist(lapply(yvals2, function(x) colMeans(x)))
ydiflocs <- unlist(lapply(yvals2, function(x) mean(x)))

pdf("/users/ndyjack/Dist_Proj/figures_final/supfig_6.1.pdf",width=8,height=12)

 #draw heatmap
  plot.new()
  par(new = "TRUE",plt = c(0.03,0.93,0.05,0.88),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:nc){
    for(j in 1:nr){
      rect(xleft=k-0.5,xright=k+0.5,ybot=yvals[1,j],ytop=yvals[2,j],border='black',col=col_pal_fxn(res[j,k]))
      text(x=k,y=mean(yvals[,j]),labels=formatC(res[j,k],digits=2,format='f'),cex=0.6,col=ifelse(res[j,k]>txtcut,'black','white'))
    }
  }
  text(x=1:nc,y=0.0,labels=colnames(res),cex=0.8)
  text(x=0.0,y=yparlocs,labels=dataset_nmz,cex=0.5,srt=0)
  text(x=nc+0.8,y=ydiflocs,srt=90,cex=0.7,labels=rep(d_vec,3))

  #color legend
  par(new = "TRUE",plt = c(0.20,0.80,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))
  mtext(side=3,text='1-G+', cex=1.2)

dev.off()

rownames(res) <- paste0("splatsim_",rownames(res))
write.csv(res,"/users/ndyjack/Dist_Proj/tables_final/suptab_8.1.csv",quote=F)
