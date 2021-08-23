#options("width"=160)
#library(Seurat)
library(PoiClaClu)
library(scran)
library(heatmap3)

#load in null data

#reinhardt cd11c+ cell experiment
#1 = IL4-positive lymph node
#2 = IL4-negative lymph node
#3 = IL-4-positive Lung
#4 = IL-4-negative Lung
counts <- as.matrix(read.table("/users/ndyjack/Dist_Proj/tables/GSE125788_10x_rawData.txt",header=T,sep="\t"))
#counts <- data.frame(do.call(cbind,counts),col.names=names(counts),row.names=rownames(counts)) 
#take the 2 most disparate populations, ie, 1 and 4
counts <- counts[,grep("1|4",colnames(counts))]		     
#remove non-expressed genes to speed up computation (maybe?? lol)
counts <- counts[-which(rowSums(counts)==0),]
#color label
labcols <-  ifelse(grepl("1",colnames(counts)),"blue","red")

#counts <- counts[,c(1:50,which(labcols=="red")[1:50])]  
#labcols <- labcols[c(1:50,which(labcols=="red")[1:50])]

#compute poisson distance
dist.pois <- PoissonDistance(t(counts))$dd
#normalize with scran
sce <- SingleCellExperiment(list(counts=counts))
sce <- computeSumFactors(sce)
sce <- normalize(sce)
counts.norm <- logcounts(sce)
#compute euclidean distance on the scran-normalized counts
dist.euc <- dist(t(counts.norm))

name.vec <- c("Euclidean Dist","log(Euc. Dist)","Poisson Dist","Log(Poi. Dist)")
dist.list <- list(dist.euc,as.dist(log(as.matrix(dist.euc)+1)),dist.pois,as.dist(log(as.matrix(dist.pois)+1)))
breakscale <- lapply(dist.list, function(x) quantile(x,c(0.01,0.99)))
breakscale <- lapply(breakscale, function(x) c(x[1],seq(x[1],x[2],length.out=255),x[2]))
hc.list <- lapply(dist.list, function(x) hclust(x,method="ward.D"))
dend.list <- lapply(hc.list, function(x) as.dendrogram(x))

plotdist.list <- lapply(dist.list, function(x) {
        dis <- as.matrix(x)
        dis[upper.tri(dis)] <- t(dis)[upper.tri(dis)]
        return(dis)
})

ncol <- 1000
col <- colorRampPalette(c("royalblue4", "ghostwhite","violetred2"))(n = ncol)
breakscale <- lapply(dist.list, function(x) {
        qs <- quantile(x,c(0.05,0.95),names=F)
        c(0,seq(qs[1],qs[2],length.out=(ncol-1)),1.1*max(x))
})

k.max <- 10
wss.list <- lapply(1:length(dist.list), function(i) {
        wss <- sapply(1:k.max, function(k) {
                clust <- cutree(hc.list[[i]],k=k)
		idx <- lapply(1:k, function(j) as.vector(which(clust==j)))
		ssq <- sapply(idx, function(j){
                        tmp <- as.matrix(dist.list[[i]])[j,j]
                        tmp <- tmp[upper.tri(tmp,diag=F)]
                        sum(tmp)
                })
        })
        wss <- sapply(wss,sum)
})
wss.scaled <- sapply(wss.list,function(x) x/max(x))
dist.cols <- c("deepskyblue3","forestgreen","red","darkmagenta")

pdf("./images/lunglymph_results_v1.pdf")
for(i in 1:length(dist.list)){
#distance heatmap (add breakscaling)
        heatmap3(plotdist.list[[i]],symm=T,Rowv=dend.list[[i]],main=name.vec[i],
                revC=T,scale='none',useRaster=T,cexRow=0.01,cexCol=0.01,
                col=col,breaks=breakscale[[i]],ColSideColors=labcols)
}
#within ssq plot
plot(x=1:k.max,y=wss.scaled[,1],type="b",pch=19,frame=F,
        xlab="# clusters K", ylab="Scaled within-cluster SSQ",col=dist.cols[1])
for(i in 2:length(dist.list)){  lines(x=1:k.max,y=wss.scaled[,i],type="b",pch=19,col=dist.cols[i]) }
legend("topright",legend=name.vec,pch=19,col=dist.cols)
dev.off()

clust_res <- lapply(hc.list,function(x) table(cutree(x,k=2),labcols))
lapply(clust_res,print)
clust_acc <- sapply(clust_res, function(x) sum(diag(x))/sum(x))
print(clust_acc)

save.image("/users/ndyjack/Dist_Proj/rdata/lunglymph.rda")


sil.list <- lapply(1:length(dist.list), function(i) {
        sil <- lapply(1:k.max, function(k) {
		clust <- cutree(hc.list[[i]],k=k)
		silhouette(x=clust,dist=dist.list[[i]])
	}) 
})

#k-mediods (can take dist mat), hierarchical,
#silouhette width (can take dist mat and arbitrary cluster assignment 
