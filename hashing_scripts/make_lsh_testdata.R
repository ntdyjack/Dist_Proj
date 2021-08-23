output_dir <- "/users/ndyjack/Dist_Proj/tables/test_datasets/"

library(heatmap3)
library(Matrix)
n_clust <- 2 #number cluters
n_cellperclust <- 500 #cells per cluster
n_features <- 10000 #10000 features
feature_means <- c(1,20) #mean expression for two genes
cell_id <- as.vector(sapply(1:n_clust, function(i) rep(i,n_cellperclust))) #celltype ID vec
#simulate expression
rna_counts <- sapply(cell_id, function(i) {
  c(rpois(n=n_features/2,lambda=feature_means[i]), rpois(n=n_features/2,lambda=feature_means[-i]))
})

#png("/users/ndyjack/Dist_Proj/images/lsh_testdata_heatmap.pdf",width=600,height=1000)
#heatmap3(rna_counts,scale="none",cexCol=0.01,cexRow=2.0,
#  useRaster=T,breaks=breakscale,col=col,method="ward.D",
#  ColSideColors=ifelse(cell_id==1,"red","blue"),ColSideLabs="",RowSideLabs="")
#dev.off()

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"lsh_testdata_2cl_full_expr.mtx"))
writeLines(as.character(cell_id),paste0(output_dir,"lsh_testdata_2cl_labs.txt"))
