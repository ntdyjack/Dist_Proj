#options("width"=160)
library(colordistance)
library(dendextend)
library(Rtsne)
library(umap)

#load in distance mats
null.poidist <- readRDS("./rdata/null_poi.rds")
null.eucdist <- readRDS("./rdata/null_euc.rds")
twoc.poidist <- readRDS("./rdata/twoc_poi.rds")
twoc.eucdist <- readRDS("./rdata/twoc_euc.rds")  

twoc.colors <- ifelse(grepl("1",colnames(as.matrix(twoc.eucdist))),"red","blue")
twoc.labs <- ifelse(grepl("1",colnames(as.matrix(twoc.eucdist))),"lymph","lung") 

dist.list <- list(null.poidist,null.eucdist,twoc.poidist,twoc.eucdist)

custom.config <- umap.defaults
custom.config$input <- "dist"

embed.list <- lapply(dist.list, function(x) {
	set.seed(15612)
	res.tsne <- Rtsne(x, is_distance=T)$Y
	res.umap <- umap(as.matrix(x),custom.config)$layout
	list(tsne=res.tsne,umap=res.umap) 
})

n <- nrow(as.matrix(null.eucdist))
name.vec <- c("Null PoiDist","Null EucDist", "2Clust PoiDist", "2Clust EucDist")
col.list <- list(rep("black",n), rep("black",n),twoc.colors,twoc.colors)

pdf("./images/compare_embeddings.pdf",width=10,height=8)

par(mfrow=c(2,2))
for(i in 1:4){
	plot(embed.list[[i]]$tsne[,1],embed.list[[i]]$tsne[,2],main=paste0("tSNE ", name.vec[i]),
		col=col.list[[i]],xlab="Dim 1", ylab = "Dim 2") 
}	
for(i in 1:4){
	plot(embed.list[[i]]$umap[,1],embed.list[[i]]$umap[,2],main=paste0("UMAP ", name.vec[i]),
		col=col.list[[i]],xlab="Dim 1", ylab = "Dim 2") 
}
dev.off()


col <- colorRampPalette(c("royalblue4", "ghostwhite","violetred2"))(n = 299)
hc.list <- lapply(dist.list, function(x) hclust(x,method="ward.D"))

dend.list <- lapply(1:4, function(i) {
	den <- as.dendrogram(hc.list[[i]])
	labels(den) <- "-"
	labels_colors(den) <- col.list[[i]][hc.list[[i]]$order]
	return(den)
}
) 
plotdist.list <- lapply(dist.list, function(x) {
	dis <- as.matrix(x)
	dis[upper.tri(dis)] <- t(dis)[upper.tri(dis)]
	return(dis)
})

pdf("./images/compare_hclust.pdf")
for(i in 1:4){
	heatmap3(plotdist.list[[i]],symm=T,Rowv=dend.list[[i]],main=name.vec[i],
		revC=T,scale='none',ColSideCol=col.list[[i]],useRaster=T,cexRow=0.01,cexCol=0.01)
}
dev.off()

lapply(hc.list[3:4], function(x) table(twoc.labs,cutree(x,2)) )

#twoc.labs   1   2
#    lung  152 678
#    lymph 674  38

#twoc.labs   1   2
#    lung   42 788
#    lymph 687  25

#read Zhang Paper
