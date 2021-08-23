#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)
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


#load in results files, 
results_dir <- "/users/ndyjack/Dist_Proj/tables/gsr/"
gsr_vals <- sapply(dists_vec, function(x) sapply(gsets, function(y) {
  tmp_nm <- paste0(results_dir,y,'.',x,'.txt')
  tmp <- as.numeric(readLines(tmp_nm))
  return(tmp)
}))




#metrics_avgs <- as.data.frame(cbind(colMeans(ari_means),colMeans(kacc_means),colMeans(gplus_means),colMeans(cor_means)))
#colnames(metrics_avgs) <- metric_nmz
#rownames(metrics_avgs) <- names(dists_vec) 
#metrics_by_gset[[5]] <- metrics_avgs
#gsets[5] <- "_avg"
#palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
#palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
#                              palette3_info$maxcolors,
#                              rownames(palette3_info)))
#cols <- sample(palette3_all, length(dists_vec))

#plotting parameters
zlim <- c(0.75,1.05)
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
nrow <- length(dists_vec)
ncol <- 4
xlim <- c(0,ncol+2)
ylim <-c(0,nrow+1)
zleglocs <- c(0.01,0.25,0.50,0.75,0.99)
zleglabs <- seq(from=zlim[1],to=zlim[2],length.out=length(zleglocs))
#as.character(pretty(zlim,length(zleglocs)))
#c("0.00","0.25","0.50","0.75","1.00")
dist_gcol <- ifelse(dists_grp2=="Probabilistic",'red','blue')
dist_pch <- ifelse(dists_grp2=="Probabilistic",15,16)

pdf("/users/ndyjack/Dist_Proj/images/gsr_heatmap.pdf",width=8,height=12)
  #draw heatmap
  tmp <- t(gsr_vals) 
  #metrics_by_gset[[i]]
  o <- order(rowMeans(tmp),decreasing=F)
  tmp <- tmp[o,]
  plot.new()
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.80),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA) 
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=nrow+0.7,labels=colnames(tmp) )
  text(x=ncol+1.0,y=(1:nrow),labels=rownames(tmp),col=dist_gcol[o])

  #color legend
  par(new = "TRUE",plt = c(0.15,0.65,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))

  #dist type legend
  par(new = "TRUE",plt = c(0.70,0.90,0.85,0.90),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=rep(0.5,2),y=c(0.2,0.8),labels=c("Probabilistic","Geometric"),col=c("red","blue"),cex=1.5)

  #title
  par(new = "TRUE",plt = c(0.20,0.80,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  text(x=0.5,y=0.5,labels="Gap Statistic Ratio",cex=2.5)

dev.off()

