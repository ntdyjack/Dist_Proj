library(PoiClaClu)
library(scran)
library(heatmap3)
library(DropletUtils)

#read in counts
counts.lunglymph <- as.matrix(read.table("/users/ndyjack/Dist_Proj/tables/GSE125788_10x_rawData.txt",header=T,sep="\t"))
counts.cd14 <-  read10xCounts("/users/ndyjack/Dist_Proj/tables/cd14_monocytes_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
counts.cd14 <- as.matrix(counts(counts.cd14))
counts.jur293t <- read10xCounts("/users/ndyjack/Dist_Proj/tables/jurkat_293t_50_50_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
counts.jur293t <- as.matrix(counts(counts.jur293t))
counts.list <- list(counts.lunglymph,counts.cd14,counts.jur293t)
counts.nz.list <- lapply(counts.list, function(x) x[-which(rowSums(x)==0),])
counts.list <- unlist(list(counts.list, counts.nz.list),recursive=F)
#counts.list <- unlist(list(counts.nz.list,counts.nz.list),recursive=F)

#alphas <-  lapply(counts.list, function(x) PoiClaClu::FindBestTransform(t(x)))
#1.0000000 1.0000000 1.0000000 0.6767347 0.6969388 0.7979592

#summary stats
ngene.vec <- sapply(counts.list, nrow)
ncell.vec <- sapply(counts.list,ncol)
N.vec <- sapply(counts.list, function(x) median(colSums(x)))
#rowmeans.list <- lapply(counts.list,function(x) log(rowMeans(x)))
#rowmeans.list[1:3] <- lapply(rowmeans.list[1:3], function(x) {
#	tmp <- x
#	tmp[is.infinite(x)] <- min(x[!is.infinite(x)])
#	tmp
#})
rowmeans.list <- unlist(list(lapply(counts.list[1:3], function(x) log(rowMeans(x)+.001)),
			lapply(counts.list[4:6], function(x) log(rowMeans(x)))),recursive=F) 
#sapply(1:length(counts.list),function(i) ifelse(i>3,log(rowMeans(counts.list[[i]])),log(rowMeans(counts.list[[i]])+1)) )
pctzero.list <- sapply(counts.list, function(x) apply(x,1, function(x) length(which(x==0))/length(x)))
pctzero.vec <- sapply(pctzero.list, function(x) sum(x==1)/length(x))
pctzero.vec <- paste0(formatC(pctzero.vec,digits=2,format="f"),"% zero")
title.vec <- c("Lung/Lymph CD4+/- (full)","CD14+ Monocytes (full)","Jurkat/293t (full)",
		"non-zero genes","non-zero genes)","non-zero genes")
subtitle.vec <- sapply(1:length(counts.list), function(i) paste(ngene.vec[i]," genes x ", ncell.vec[i]," cells"))
yticks <- seq(0,1,by=0.20)
ylabs <- formatC(yticks,digits=1,format="f")
xticks.list <- lapply(rowmeans.list, function(x) seq(min(x),max(x),length.out=5))
#xticks.list <- unlist(list(xticks.list,xticks.list),recursive=F)
xlabs.list <- lapply(xticks.list, function(x) formatC(x,digits=2,format="f"))

#accessory functions
predict_zeros_poi<-function(x){ exp(-exp(x)) }
predict_zeros_binom<-function(x,N){ (1-exp(x)/N)^N }
predict_zeros_nb<-function(x,phi=2){ exp(-phi*log1p(exp(x-log(phi)))) }

#plotting
png("./images/frac_zero_plots.png",width=1200,height=800)
par(mfrow=c(2,3),mar=c(4.1,5.1,4.1,1.1))
for(i in 1:length(counts.list)){
	xlab <- ifelse(i>3,"log(mean expression) [gene]","log(mean expression+.001) [gene]")
	legloc <- ifelse(i>3,"bottomleft","topright")
	plot(x=rowmeans.list[[i]],y=pctzero.list[[i]],pch=1,cex=0.2,
		ylab="",yaxt="n",xlab="",xaxt="n",main="")
	axis(side=1,labels=xlabs.list[[i]],at=xticks.list[[i]],cex=0.8,mgp=c(3, .5, 0))
	axis(side=2,labels=ylabs,at=yticks,"fraction zero [gene]",cex=0.8,las=2,mgp=c(3, .5, 0))
	mtext(side=1,text=xlab,line=2.0,cex=1.5)
	mtext(side=2,text="fraction zero [gene]",line=2.0,cex=1.5)
	mtext(side=3,text=title.vec[i],line=2.0,cex=2.0)
        mtext(side=3,text=subtitle.vec[i],line=0.5,cex=1.0)
	text(x=xlabs.list[[i]][2],y=yticks[5],labels=pctzero.vec[i],cex=2)

	xlo <- min(rowmeans.list[[i]])
	xhi <- max(rowmeans.list[[i]])
	curve(predict_zeros_binom(x,N=N.vec[i]),from=xlo,to=xhi,col="blue",lwd=2,add=TRUE)
	curve(predict_zeros_poi,from=xlo,to=xhi,col="green",lwd=2,lty=2,add=TRUE)
	curve(predict_zeros_nb(x,phi=4),from=xlo,to=xhi,col="red",lwd=3,lty=3,add=TRUE)
	legend(legloc,c("Multinomial","Poisson","Negative Binomial"),
		lty=c(1,2,3),lwd=c(4,3,3),col=c("blue","green","red"))
}
dev.off()

