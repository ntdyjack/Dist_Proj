####################
#Custum Violin Plot#
####################
library(beeswarm)
ggdat_pch <- c(16,16,16,15,15)
atl <- 1:5
xlims <- c(0.5,5.5)
for(i in 1:l){
    j <- as.character(i-1)
    degs <- sort(linear_test_byclust[[i]]$gene[linear_test_byclust[[i]]$fdr<5e-100])
    print(paste0("plotting... cluster",j,", n=",length(degs)," DEGs"))
    #don't make a plot if there's no degs that meet bonferroni
    if(length(degs)==0){next}
    tmp <- SubsetData(raw,ident.use=j)
    tmp_id <- tmp$final
    tmp <- as.matrix(tmp$RNA@data)
    pdf(paste0("./IMAGES/combinedFinal_cluster",j,"_LMgenes_5e100.pdf"),width=10,height=8,family='sans')
    par(mar=c(1.1,3.1,1.6,0.1))
    for(x in degs){
        print(x)
        ggdat <- lapply(final_ref, function(y) tmp[x,tmp_id==y])
        ggident_use <- which(sapply(ggdat,length)>1 & !sapply(ggdat,function(x) all(x==x[1])))
        ggdens <- lapply(1:length(ggdat), function(i) NA)
        ggdens[ggident_use] <- lapply(ggdat[ggident_use], function(y) density(y,na.rm=F,from=min(y),to=max(y),bw=0.3))
        ylims <- c(round(min(unlist(ggdat)),2),round(max(unlist(ggdat)),2))
        ytxl <- seq(ylims[1],ylims[2],length.out=5)
        plot(x=0,y=0,type='n',xaxt='n',yaxt='n',xlim=xlims,ylim=ylims,xlab="",ylab="")
        for(k in ggident_use){
            freqs <- ggdens[[k]]$y
            freqs <- freqs/(max(freqs)*2.5)
            exprs <- ggdens[[k]]$x
            exprs <- c(exprs,rev(exprs))
            freqs <- c(k+freqs,rev(k-freqs))
            polygon(x=freqs,y=exprs,col=final_color_ref[k],lwd=2,border=NA)
        }
        beeswarm(ggdat,pch=ggdat_pch,at=atl,add=T,cex=1.0,corral="omit",corralWidth=0.5)
        mtext(text=paste0(x,", FDR=",formatC(linear_test_byclust[[i]][x,'fdr'],digits=2,format="e")),side=3,cex=2.0,line=0.0)
        axis(side=2,labels=formatC(ytxl,digits=1,format="f"),cex=0.7,line=0,at=ytxl,las=2,mgp=c(0.0,0.6,0))
        mtext(text="Scaled Expression",side=2,cex=1.5,line=1.8)
        mtext(text=final_ref,side=1,cex=1.5,line=0.1,at=atl,mgp=c(0.0,0.4,0))
        pctzero <- formatC(sapply(ggdat, function(x) length(which(x==0))/length(x))*100,digits=0,format="f")
        text(x=atl,y=ifelse(ylims[1]==0,ylims[1]-(0.02*ylims[2]),0.98*ylims[1]),labels=paste0(pctzero,"% Zero"),cex=0.7)
    }
    dev.off()
}

