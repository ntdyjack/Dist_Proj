#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
library(pacman)
p_load(heatmap3,beeswarm,RColorBrewer,viridis,leaflet)

map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

#trajcor
files = c(
  "cellmix_rnaseq_traj",
  "rnamixcs2_rnaseq_traj",
  "rnamixss_rnaseq_traj"
)

gsets = c(
  "_hvg1k"
)

datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

#names of the distances being tested
dists_vec <- c(
  'L1'=0, 'L2'=1, 'SQL2'=2,
  'JSM'=3, 'JSD'=4, 'MKL'=5, 
  'HEL'=6, 'BCM'=7, 'BCD'=8,
  'TVD'=9, 'LLR'=10, 'UWLLR'=12,
  'ISD'=13, 'RISD'=14, 'SIS'=23
)
l <- length(dists_vec)

#trajcor real data
results_dir <- '/users/ndyjack/Dist_Proj/tables/trajcor_new/'
cor_sorted <- results_sorted <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- read.table(paste0(results_dir,y,".",x,".txt"),header=T)
  return(mean(tmp$pearson))
}))
cor_rna_final <- colMeans(cor_sorted)

#trajcor simulated data
res_dir <- "/users/ndyjack/Dist_Proj/tables/trajcor_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_tc <- sapply(res_files, function(x) {
  tmp <- read.table(paste0(res_dir,x),header=T)
  tmp <- mean(tmp$pearson)
  return(tmp)
})
cor_sorted <- sapply(dists_vec, function(x) {
  tmp <- res_tc[grep(paste0('\\.',x,'\\.txt'),names(res_tc))]
})
cor_sim_final <- colMeans(cor_sorted)


# real data
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
datasets_vec <- unlist(lapply(files, function(x) paste0(x,gsets)))

#ARI real data
results_dir <- "/users/ndyjack/Dist_Proj/tables/ari_new/"
ari_vals <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".kmeds.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
ari_rna_final <- colMeans(ari_vals)

#ARI sim data
res_dir <- "/users/ndyjack/Dist_Proj/tables/ari_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
ari_sim_final <- sapply(dists_vec, function(x) {
  tmpn <- grep(paste0('\\.',x,'\\.kmeds\\.txt'),res_files,value=T)
  tmpv <- sapply(tmpn, function(y) {
    tmps <- readLines(paste0(res_dir,y))
    tmps <- as.numeric(strsplit(tmps,' ')[[1]][2])
    return(tmps)
  })
  return(mean(tmpv))
})

#kAcc real data
results_dir <- "/users/ndyjack/Dist_Proj/tables/kacc_new/"
kacc_vals <- sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".k30.txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
kacc_rna_final <- colMeans(kacc_vals)

#kAcc sim data
res_dir <- "/users/ndyjack/Dist_Proj/tables/kacc_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
res_files <- grep('k30', res_files, value=T)
kacc_sim_final <- sapply(dists_vec, function(x) {
  tmpn <- grep(paste0('\\.',x,'.k30.txt'),res_files,value=T)
  tmpv <- sapply(tmpn, function(y) {
    tmps <- readLines(paste0(res_dir,y))
    tmps <- as.numeric(strsplit(tmps,' ')[[1]][2])
    return(tmps)
  })
  return(mean(tmpv))
})

#1-G+ real data
results_dir <- "/users/ndyjack/Dist_Proj/tables/gplus_new/"
gplus_sorted <- 1 - sapply(dists_vec, function(x) sapply(datasets_vec, function(y) {
  tmp <- readLines(paste0(results_dir,y,".",x,".txt"))
  tmp <- as.numeric(strsplit(tmp," ")[[1]][2])
}))
gplus_rna_final <- colMeans(gplus_sorted)


#1-G+ sim data
res_dir <- "/users/ndyjack/Dist_Proj/tables/gplus_new/"
res_files <- list.files(res_dir)
res_files <- grep('splat',res_files,value=T)
gplus_sim_final <- 1 - sapply(dists_vec, function(x) {
  tmpn <- grep(paste0('\\.',x,'.txt'),res_files,value=T)
  tmpv <- sapply(tmpn, function(y) {
    tmps <- readLines(paste0(res_dir,y))
    tmps <- as.numeric(strsplit(tmps,' ')[[1]][2])
    return(tmps)
  })
  return(mean(tmpv))
})


#GS real data
results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/hvg1k.'
gs_rna_final <- sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  tmp <- map_fxn(tmp,0,1,min(tmp),max(tmp))
  return(tmp[1]-tmp[2])
})

#GS sim data
results_dir <- '/users/ndyjack/Dist_Proj/tables/gapstat/splat.hvg1k.'
gs_sim_final <- sapply(dists_vec, function(y) {
  tmp <- paste0(results_dir,y,'.txt')
  tmp <- readLines(tmp)
  tmp <- as.numeric(strsplit(tmp," ",)[[1]])
  tmp <- map_fxn(tmp,0,1,min(tmp),max(tmp))
  return(tmp[1]-tmp[2])
})


sim_plot <- list(gs_sim_final,kacc_sim_final, gplus_sim_final, ari_sim_final, cor_sim_final)
sim_plot <- t(sapply(sim_plot, function(x) map_fxn(x,0,1,min(x),max(x))))
sim_fin <- colMeans(sim_plot)
o_sim <- order(sim_fin)

rna_plot <- list(gs_rna_final,kacc_rna_final, gplus_rna_final, ari_rna_final, cor_rna_final)
rna_plot <- t(sapply(rna_plot, function(x) map_fxn(x,0,1,min(x),max(x))))
rna_fin <- colMeans(rna_plot)
o_rna <- order(rna_fin)


rnk_plot <- rbind(rna_fin,sim_fin)[,o_rna]
rnk_plot <- t(apply(rnk_plot,1,function(x) map_fxn(x,0,1,min(x),max(x))))

#rna_plot <- do.call(rbind( lapply(list(gs_rna_final,kacc_rna_final, gplus_rna_final, ari_rna_final, trajcor_rna_final) ),  function(x) map_fxn(tmp,0,1,min(x),max(x))))
#rna_fin <- rowMeans(sim_plot)

#plotting parameters
palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
palette3_all <- unlist(mapply(brewer.pal,palette3_info$maxcolors,rownames(palette3_info)))
set.seed(1254)
col_labs <- sample(palette3_all, length(dists_vec))
nlevs <- l
zlim <- c(0,1)
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
#xlim <- c(0,l+2)
#ylim <-c(0,1)
pch_labs <- rep(c(15,16,17),5)
zleglocs <- c(0.01,0.99)
zleglabs <- c("Worst","Best")
stat_nmz <- c('GapStat','kAcc','1-G+','ARI','TrajCor')
ncol <- l
nrow <- length(stat_nmz)
xlim <- c(0,ncol+1)
ylim <- c(0,nrow+1)

pdf("/users/ndyjack/Dist_Proj/figures_final/figure_5.1.pdf",width=8,height=12)
  plot.new() 
 
  #color bar legend for library size
  par(new = "TRUE",plt = c(0.20,0.80,0.90,0.97),las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  mtext(side=3,text='Mean Performance', cex=1.2)
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.4,las=1,mgp=c(1.0, .4, 0))


  #GS(r?) real data (?)
  #1-G+ real data
  #ARI real data
  #kAcc real data
  #trajcor real data
  #heatmap (bottom)
  tmp <- rna_plot[,o_rna] 
  par(new = "TRUE",plt = c(0.03,0.97,0.55,0.85),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=0.3,labels=names(dists_vec)[o_rna],cex=0.7)
  text(x=0.2,y=(1:nrow),labels=stat_nmz,srt=90)
  text(x=xlim[2]-0.2,y=mean(ylim),labels='scRNA-seq (hvg1k)',cex=1.0,srt=90)

  #GS(r?) sim data (?)
  #1-G+ sim data
  #ARI sim data
  #kAcc sim data
  #trajcor real data  
  tmp <- sim_plot[,o_sim]
  par(new = "TRUE",plt = c(0.03,0.97,0.25,0.55),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=(1:ncol),y=0.3,labels=names(dists_vec)[o_sim],cex=0.7)
  text(x=0.2,y=(1:nrow),labels=stat_nmz,srt=90)
  text(x=xlim[2]-0.2,y=mean(ylim),labels='Splatter (hvg1k)',srt=90)

  #overall ranking (?)
  nrow <- 2
  ylim <- c(0,nrow+1)
  tmp <- rnk_plot
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.25),las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(j in 1:nrow){
    for(k in 1:ncol){
      rect(xleft=k-0.5,xright=k+0.5,ybot=j-0.5,ytop=j+0.5,border='black',col=col_pal_fxn(tmp[j,k]))
    }
  }
  text(x=seq(1,l,by=2),y=0.15,labels=names(dists_vec)[o_rna][seq(1,l,by=2)],cex=1.2)
  text(x=seq(2,l-1,by=2),y=0.35,labels=names(dists_vec)[o_rna][seq(2,l-1,by=2)],cex=1.2)
  text(x=xlim[2]-0.2,y=(1:nrow),labels=rev(c('SIM','RNA')),srt=90)
  text(x=mean(xlim),y=ylim[2]-0.3,labels='Rank',cex=1.2)

dev.off()


