#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet,RcppCNPy,reticulate,Matrix)
np <-import("numpy")
rndup <- function(x,d) { ceiling(x*d) / d }
rnddn <- function(x,d) { floor(x*d) / d }


#####
gsets = c(
  "full",
  "hvg5k",
  "hvg1k",
  "hvg500"
)

#names of the distances being tested
dists_vec <- c(
  'L2'=1,
  'JSM'=3, 
  'LLR'=10
)

#calculate library size factors
input_dir  <-  "/users/ndyjack/Dist_Proj/tables/test_datasets/only293t_rnaseq_1cl_"
ls_vals <- lapply(gsets, function(x) {
  nm <- paste0(input_dir,x,'_expr.mtx')
  count_mat <- readMM(nm)
  cs_vec <- colSums(count_mat)
  return(cs_vec)
})

#read in PAM labels
input_dir  <-  "/users/ndyjack/Dist_Proj/tables/test_datasets/only293t_rnaseq_1cl_"
pam_vals <- lapply(gsets, function(x) lapply(dists_vec, function(y) {
  nm <- paste0(input_dir,x,'.',y,'.labs.txt')
  labs <- readLines(nm)
  cols <- ifelse(labs==1,'red','blue')
  #count_mat <- readMM(nm)
  #cs_vec <- colSums(count_mat)
  return(cols)
}))

#load in results files, 
results_dir <- "/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_"
dr_vals <- lapply(gsets, function(x) lapply(dists_vec, function(y) {
  nm <- paste0(results_dir,x,'.',y,'.npy')
  dist_mat <- np$load(nm)
  dist_mat <- as.dist(dist_mat)
  mds_mat <- as.matrix(cmdscale(d=dist_mat,k = 2))
  #tmp <- as.numeric(readLines(tmp_nm))
  return(mds_mat)
}))

plot_lims <- list(
  c(0.20,0.95,0.10,0.30),
  c(0.20,0.95,0.40,0.60),
  c(0.20,0.95,0.70,0.90)
)

pdf("/users/ndyjack/Dist_Proj/images/null_dr_mds_libsize.pdf",width=4,height=12)
for(i in 1:length(gsets)){

  #title
  plot.new()
  par(new = "TRUE",plt = c(0.20,0.80,0.92,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.0)

  #point colors
  zlim <- c(min(ls_vals[[i]]),max(ls_vals[[i]]))
  nlevs <- 100
  levs <- seq(0,1,length.out=nlevs)
  col_pal <- heat.colors(nlevs-1) #viridis(nlevs-1)
  col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
  col_vec <- col_pal_fxn(ls_vals[[i]])
  col_vec <- paste0(col_vec,'64')


  for(j in 1:length(dists_vec)){
  #scatter plot
    par(new = "TRUE",plt = plot_lims[[j]],las = 1,cex.axis = 1)
    tmp <- dr_vals[[i]][[j]]
    xval <- tmp[,1]
    yval <- tmp[,2]
    xlim <- c(min(xval),max(xval))
    ylim <- c(min(yval),max(yval))
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
    box()
    points(x=xval,y=yval,col=col_vec,cex=0.8,pch=16)
    mtext(side=1,text='MDS 1',cex=1.2,line=1.5,las=1)
    mtext(side=2,text='MDS 2',cex=1.2,line=2.5,las=3)
    mtext(side=3,text=names(dists_vec)[j],cex=2.0,line=0.3)
  }
}
dev.off()

#color bar
#pdf(""

pdf("/users/ndyjack/Dist_Proj/images/null_dr_mds_pamclust.pdf",width=4,height=12)
for(i in 1:length(gsets)){

  #title
  plot.new()
  par(new = "TRUE",plt = c(0.20,0.80,0.92,0.99),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels=sub("_","",gsets[i]),cex=3.0)

  for(j in 1:length(dists_vec)){
  #scatter plot
    par(new = "TRUE",plt = plot_lims[[j]],las = 1,cex.axis = 1)
    tmp <- dr_vals[[i]][[j]]
    xval <- tmp[,1]
    yval <- tmp[,2]
    xlim <- c(min(xval),max(xval))
    ylim <- c(min(yval),max(yval))
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i",yaxs = "i")
    box()
    points(x=xval,y=yval,col=pam_vals[[i]][[j]],cex=0.8,pch=16)
    mtext(side=1,text='MDS 1',cex=1.2,line=1.5,las=1)
    mtext(side=2,text='MDS 2',cex=1.2,line=2.5,las=3)
    mtext(side=3,text=names(dists_vec)[j],cex=2.0,line=0.3)
  }
}
dev.off()
