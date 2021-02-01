#ml old_conda_R/3.5.x
options("width"=180)
library(pacman)
p_load(DropletUtils,heatmap3,viridis,readxl,Matrix,scran)

ncol <- 1000
col <- viridis(ncol)

output_dir <- "/users/ndyjack/Dist_Proj/tables/test_datasets/"

#############################
#######cite seq 10k pbmc data
#############################
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/pbmc_10k_citeseq") 
rna_counts <- as.matrix(counts(sce))
rna_colsums <- colSums(rna_counts[grep("ENSG",rownames(rna_counts)),],)
#mito genes
mito_genes <- rownames(rna_counts)[grep("^MT-",sce@rowRanges@elementMetadata@listData$Symbol)]
pct_mito <- colSums(rna_counts[mito_genes,]) / rna_colsums
ptn_counts <- rna_counts[!grepl("ENSG",rownames(rna_counts)),]
ptn_colsums <- colSums(ptn_counts)
rna_counts <- rna_counts[grepl("ENSG",rownames(rna_counts)),]

#qc data
qtls_dropm <- quantile(pct_mito,c(0.005,0.995))
qtls_dropr <- quantile(rna_colsums,c(0.005,.995))
qtls_dropp <- quantile(ptn_colsums,c(0.005,.995))
cells_drop <- which(
  ptn_colsums > qtls_dropp[2] |  ptn_colsums < qtls_dropp[1] |
  rna_colsums > qtls_dropr[2] |  rna_colsums < qtls_dropr[1] |
  pct_mito > qtls_dropm[2] |  pct_mito < qtls_dropm[1]
)

pdf("/users/ndyjack/Dist_Proj/images/pbmc_citeseq_qc.pdf",width=12,height=8)
par(mfrow=c(2,2))
for(x in rownames(ptn_counts)){
  hist(breaks=100,ptn_counts[x,],main=x)
}
hist(breaks=100,ptn_colsums)
abline(v=qtls_dropp,col="red",lty=2)
hist(breaks=100,rna_colsums)
abline(v=qtls_dropr,col="red",lty=2)
hist(breaks=100,pct_mito)
abline(v=qtls_dropm,col="red",lty=2)
dev.off()

#subset to cells that pass QC, normalize
rna_counts <- rna_counts[,-cells_drop]
ptn_counts <- ptn_counts[,-cells_drop]
ptn_unit <- apply(ptn_counts,2,function(x) x/sum(x))
ptn_norm <- apply(ptn_counts,2,function(x) log( ((x/sum(x)) + 1)))

prots_use <- c("CD3","CD16","CD8a","CD19","CD4","CD14")
dat_plot <- ptn_norm[prots_use,]
dat_plot <- apply(dat_plot, 1, function(x) (x - mean(x)) / sd(x))
qtls <- quantile(dat_plot,c(0.05,0.95))
breakscale <- c(min(dat_plot),seq(qtls[1],qtls[2],length.out=ncol-1),max(dat_plot))

hc <- hclust(dist(dat_plot),method="ward.D")
clust <- cutree(hc,5)
clust_colref <- c("red","blue","darkgreen","darkorange","darkturquoise")
clust_colvals <- sapply(clust, function(i) clust_colref[i])
dend <- as.dendrogram(hc)
table(clust_colvals)

pdf("/users/ndyjack/Dist_Proj/images/pbmc_citeseq_heatmap.pdf",width=10,height=8)
heatmap3(t(dat_plot),scale="none",cexCol=0.01,cexRow=2.0,
  useRaster=F,breaks=breakscale,col=col,method="ward.D",
  Colv = dend, ColSideColors=clust_colvals,ColSideLabs="")
dev.off()

#compute variant genes
rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),]
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"pbmc_citeseq_5cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"pbmc_citeseq_5cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"pbmc_citeseq_5cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"pbmc_citeseq_5cl_hvg5k_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"pbmc_citeseq_5cl_labs.txt"))


cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

##########################
####cite seq 10k malt data
##########################
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/malt_10k_citeseq")
rna_counts <- as.matrix(counts(sce))
rna_colsums <- colSums(rna_counts[grep("ENSG",rownames(rna_counts)),],)
pct_mito <- colSums(rna_counts[mito_genes,]) / rna_colsums
ptn_counts <- rna_counts[!grepl("ENSG",rownames(rna_counts)),]
ptn_colsums <- colSums(ptn_counts)
rna_counts <- rna_counts[grepl("ENSG",rownames(rna_counts)),]

#qc data
qtls_dropm <- quantile(pct_mito,c(0.005,0.995))
qtls_dropr <- quantile(rna_colsums,c(0.005,.995))
qtls_dropp <- quantile(ptn_colsums,c(0.005,.995))
cells_drop <- which(
  ptn_colsums > qtls_dropp[2] |  ptn_colsums < qtls_dropp[1] |
  rna_colsums > qtls_dropr[2] |  rna_colsums < qtls_dropr[1] |
  pct_mito > qtls_dropm[2] |  pct_mito < qtls_dropm[1]
)

pdf("/users/ndyjack/Dist_Proj/images/malt_citeseq_qc.pdf",width=12,height=8)
par(mfrow=c(2,2))
for(x in rownames(ptn_counts)){
  hist(breaks=100,ptn_counts[x,],main=x)
}
hist(breaks=100,ptn_colsums)
abline(v=qtls_dropp,col="red",lty=2)
hist(breaks=100,rna_colsums)
abline(v=qtls_dropr,col="red",lty=2)
hist(breaks=100,pct_mito)
abline(v=qtls_dropm,col="red",lty=2)
dev.off()

#subset to keep cells, normalize
rna_counts <- rna_counts[,-cells_drop]
ptn_counts <- ptn_counts[,-cells_drop]
ptn_unit <- apply(ptn_counts,2,function(x) x/sum(x))
ptn_norm <- apply(ptn_counts,2,function(x) log( ((x/sum(x)) + 1)))

prots_use <- c("CD3","CD19") 
dat_plot <- ptn_norm[prots_use,]
dat_plot <- apply(dat_plot, 1, function(x) (x - mean(x)) / sd(x))
qtls <- quantile(dat_plot,c(0.05,0.95))
breakscale <- c(min(dat_plot),seq(qtls[1],qtls[2],length.out=ncol-1),max(dat_plot))

hc <- hclust(dist(dat_plot),method="ward.D")
clust <- cutree(hc,2)
clust_colref <- c("red","blue")
clust_colvals <- sapply(clust, function(i) clust_colref[i])
dend <- as.dendrogram(hc)
table(clust_colvals)

pdf("/users/ndyjack/Dist_Proj/images/malt_citeseq_heatmap.pdf",width=10,height=8)
heatmap3(t(dat_plot),scale="none",cexCol=0.01,cexRow=2.0,
  useRaster=F,breaks=breakscale,col=col,method="ward.D",
  Colv = dend, ColSideColors=clust_colvals,ColSideLabs="")
dev.off()


#compute variant genes
rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),]
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"malt_citeseq_2cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"malt_citeseq_2cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"malt_citeseq_2cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"malt_citeseq_2cl_hvg5k_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"malt_citeseq_2cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

##########################
####jurkat 293t 50/50 data
##########################

#load in data from jurkat 293t 50/50
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/jurkat_293t_50_50_filtered_gene_bc_matrices/filtered_matrices_mex/hg19")
#make identity vector for jurkat vs 239t
clust <- read_excel("/users/ndyjack/Dist_Proj/tables/41467_2017_BFncomms14049_MOESM830_ESM.xlsx")
clust <- as.integer(as.factor(clust$`Species Call by Expression`))
rna_counts <- as.matrix(counts(sce))
rna_colsums <- colSums(rna_counts[grep("ENSG",rownames(rna_counts)),],)
pct_mito <- colSums(rna_counts[mito_genes,]) / rna_colsums
#remove zero-expression cite-seq cell
qtls_dropm <- quantile(pct_mito,c(0.005,0.995))
qtls_dropr <- quantile(rna_colsums,c(0.005,.995))
cells_drop <- which(
  rna_colsums > qtls_dropr[2] |  rna_colsums < qtls_dropr[1] |
  pct_mito > qtls_dropm[2] |  pct_mito < qtls_dropm[1]
)

pdf("/users/ndyjack/Dist_Proj/images/jurkat293t_qc.pdf",width=12,height=8)
hist(breaks=100,rna_colsums)
abline(v=qtls_dropr,col="red",lty=2)
hist(breaks=100,pct_mito)
abline(v=qtls_dropm,col="red",lty=2)
dev.off()


#compute variant genes
rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),]
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"jurkat293t_rnaseq_2cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"jurkat293t_rnaseq_2cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"jurkat293t_rnaseq_2cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"jurkat293t_rnaseq_2cl_hvg5k_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"jurkat293t_rnaseq_2cl_labs.txt"))


cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

##########################
###CellBench datasets
##########################
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767
# wget https://github.com/Shians/scBenchData/, https://github.com/LuyiTian/sc_mixology/tree/master/data/csv
#tmp <- load("/users/ndyjack/Dist_Proj/tables/CellBench/sincell_with_class.RData")
#library(BiocManager)
#BiocManager::install("CellBench")
#library(CellBench)
#datasets <- CellBench::load_sc_data()
#CellBench::load_cell_mix_data,CellBench::load_mrna_mix_data,load_sc_data()

#10x cellbench 5cluster data
cellbench_dir <- "/users/ndyjack/Dist_Proj/tables/CellBench/LuyiTian/data/csv"
#cellbench_files <- list.files(cellbench_dir)
#rna_counts_geo <- read.csv("/users/ndyjack/Dist_Proj/tables/CellBench/GEO/GSM3618014_gene_count.csv",row.names=1)
meta_data <- read.csv(paste0(cellbench_dir,"/sc_10x_5cl.metadata.csv"))
clust <- as.integer(as.factor(meta_data$cell_line))
rna_counts <- as.matrix(read.csv(paste0(cellbench_dir,"/sc_10x_5cl.count.csv")))
#> dim(rna_counts_git)
#[1] 11786  3918
#> dim(rna_counts_geo)
#[1] 32895  5001
all.equal(rownames(meta_data),colnames(rna_counts))

#compute variant genes
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"cellbench10x_rnaseq_5cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,file="cellbench10x_rnaseq_5cl_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,file="cellbench10x_rnaseq_5cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"/cellbench10x_rnaseq_5cl_hvg500_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"cellbench10x_rnaseq_5cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#10x cellbench 3cluster data
meta_data <- read.csv(paste0(cellbench_dir,"/sc_10x.metadata.csv"))
clust <- as.integer(as.factor(meta_data$cell_line))
rna_counts <- as.matrix(read.csv(paste0(cellbench_dir,"/sc_10x.count.csv")))
all.equal(rownames(meta_data),colnames(rna_counts))

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"cellbench10x_rnaseq_3cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"cellbench10x_rnaseq_3cl_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"cellbench10x_rnaseq_3cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"cellbench10x_rnaseq_3cl_hvg500_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"cellbench10x_rnaseq_3cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#cllbench dropseq 3cluster data
meta_data <- read.csv(paste0(cellbench_dir,"/sc_dropseq.metadata.csv"))
clust <- as.integer(as.factor(meta_data$cell_line))
rna_counts <-as.matrix(read.csv(paste0(cellbench_dir,"/sc_dropseq.count.csv")))
all.equal(rownames(meta_data),colnames(rna_counts))

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"cellbenchDS_rnaseq_3cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"cellbenchDS_rnaseq_3cl_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"cellbenchDS_rnaseq_3cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"cellbenchDS_rnaseq_3cl_hvg500_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"cellbenchDS_rnaseq_3cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#cllbench cellseq2 3cluster data
meta_data <- read.csv(paste0(cellbench_dir,"/sc_celseq2.metadata.csv"))
clust <- as.integer(as.factor(meta_data$cell_line))
rna_counts <-as.matrix(read.csv(paste0(cellbench_dir,"/sc_celseq2.count.csv")))
all.equal(rownames(meta_data),colnames(rna_counts))
rna_counts<- rna_counts[-which(rowSums(rna_counts)==0),]

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"cellbenchCS2_rnaseq_3cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"cellbenchCS2_rnaseq_3cl_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"cellbenchCS2_rnaseq_3cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"cellbenchCS2_rnaseq_3cl_hvg500_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"cellbenchCS2_rnaseq_3cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#cellbench cellseq2 5cluster data
rna_list <- lapply(c("/sc_celseq2_5cl_p1.count.csv","/sc_celseq2_5cl_p2.count.csv","/sc_celseq2_5cl_p3.count.csv"), function(x) read.csv(paste0(cellbench_dir,x)))
common_genes <- unique(unlist(lapply(rna_list,rownames)))
common_genes <- common_genes[common_genes %in% rownames(rna_list[[1]]) & common_genes %in% rownames(rna_list[[2]]) & common_genes %in% rownames(rna_list[[3]])]
rna_list <- lapply(rna_list, function(x) x[common_genes,])
rna_counts <- as.matrix(do.call(cbind,rna_list))

meta_list <- lapply(c("/sc_celseq2_5cl_p1.metadata.csv","/sc_celseq2_5cl_p2.metadata.csv","/sc_celseq2_5cl_p3.metadata.csv"), function(x) read.csv(paste0(cellbench_dir,x)))
meta_data <- do.call(rbind,meta_list)
clust <- as.integer(as.factor(meta_data$cell_line)) #unlist(lapply(meta_list, function(x) x$'cell_line'))
all.equal(rownames(meta_data),colnames(rna_counts))
#there are some errors with the strings, but this seems OK

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"cellbenchCS2_rnaseq_5cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"cellbenchCS2_rnaseq_5cl_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"cellbenchCS2_rnaseq_5cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"cellbenchCS2_rnaseq_5cl_hvg500_expr.mtx"))
writeLines(as.character(clust),paste0(output_dir,"cellbenchCS2_rnaseq_5cl_labs.txt"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#cellbench cellmix datasets (combination of 4 plates, the 5th is just controls)
rna_list <- lapply(paste0("/",paste0("cellmix",1:4),".count.csv"), function(x) read.csv(paste0(cellbench_dir,x)))
common_genes <- unique(unlist(lapply(rna_list,rownames)))
common_genes <- common_genes[which(apply(sapply(rna_list, function(x) common_genes %in% rownames(x)),1,all))]
rna_list <- lapply(rna_list, function(x) x[common_genes,])
rna_counts <- as.matrix(do.call(cbind,rna_list))

meta_list <- lapply(paste0("/",paste0("cellmix",1:4),".metadata.csv"), function(x) read.csv(paste0(cellbench_dir,x)))
meta_list <- lapply(meta_list, function(x) x[,c('H2228','H1975','HCC827')])
meta_data <- do.call(rbind,meta_list)
all.equal(rownames(meta_data),colnames(rna_counts))
#there are some errors with the strings, but this seems OK

#remove the libraries that are only 3 cells
cells_drop  <- which(rowSums(meta_data)==3)
rna_counts <- rna_counts[,-cells_drop]
meta_data <- meta_data[-cells_drop,]

#see https://github.com/LuyiTian/sc_mixology/blob/master/script/trajectory/trajectory_RNAmix.R
group <- apply(meta_data,1,function(x) paste0(x,collapse=" "))

H2228_to_H1975 <- H2228_to_HCC827 <- H1975_to_HCC827 <- rep(NA,length(group))

H2228_to_H1975[group=="9 0 0"] = 0
H2228_to_H1975[group=="8 1 0"] = 1
H2228_to_H1975[group=="7 1 1"] = 2
H2228_to_H1975[group=="7 2 0"] = 3
H2228_to_H1975[group=="6 3 0"] = 4
H2228_to_H1975[group=="5 2 2"] = 5
H2228_to_H1975[group=="5 4 0"] = 6
H2228_to_H1975[group=="3 3 3"] = 7
H2228_to_H1975[group=="4 5 0"] = 8
H2228_to_H1975[group=="2 5 2"] = 9
H2228_to_H1975[group=="3 6 0"] = 10
H2228_to_H1975[group=="2 7 0"] = 11
H2228_to_H1975[group=="1 7 1"] = 12
H2228_to_H1975[group=="1 8 0"] = 13
H2228_to_H1975[group=="0 9 0"] = 14

H2228_to_HCC827[group=="9 0 0"] = 0
H2228_to_HCC827[group=="8 0 1"] = 1
H2228_to_HCC827[group=="7 1 1"] = 2
H2228_to_HCC827[group=="7 0 2"] = 3
H2228_to_HCC827[group=="6 0 3"] = 4
H2228_to_HCC827[group=="5 2 2"] = 5
H2228_to_HCC827[group=="5 0 4"] = 6
H2228_to_HCC827[group=="3 3 3"] = 7
H2228_to_HCC827[group=="4 0 5"] = 8
H2228_to_HCC827[group=="2 2 5"] = 9
H2228_to_HCC827[group=="3 0 6"] = 10
H2228_to_HCC827[group=="2 0 7"] = 11
H2228_to_HCC827[group=="1 1 7"] = 12
H2228_to_HCC827[group=="1 0 8"] = 13
H2228_to_HCC827[group=="0 0 9"] = 14

H1975_to_HCC827[group=="0 9 0"] = 0
H1975_to_HCC827[group=="0 8 1"] = 1
H1975_to_HCC827[group=="1 7 1"] = 2
H1975_to_HCC827[group=="0 7 2"] = 3
H1975_to_HCC827[group=="0 6 3"] = 4
H1975_to_HCC827[group=="2 5 2"] = 5
H1975_to_HCC827[group=="0 5 4"] = 6
H1975_to_HCC827[group=="3 3 3"] = 7
H1975_to_HCC827[group=="0 4 5"] = 8
H1975_to_HCC827[group=="2 2 5"] = 9
H1975_to_HCC827[group=="0 3 6"] = 10
H1975_to_HCC827[group=="0 2 7"] = 11
H1975_to_HCC827[group=="1 1 7"] = 12
H1975_to_HCC827[group=="0 1 8"] = 13
H1975_to_HCC827[group=="0 0 9"] = 14

traj_orders <- data.frame(H2228_to_H1975,H2228_to_HCC827,H1975_to_HCC827)
#meta_data[which(apply(traj_orders,1,function(x) all(is.na(x)))),]

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"cellmix_rnaseq_traj_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"cellmix_rnaseq_traj_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"cellmix_rnaseq_traj_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"cellmix_rnaseq_traj_hvg500_expr.mtx"))
write.table(traj_orders,paste0(output_dir,"cellmix_rnaseq_traj_labs.txt"),quote=F,row.names=F)

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))


#cellbench celseq2 RNA-mix 
rna_counts <-as.matrix(read.csv(paste0(cellbench_dir,"/RNAmix_celseq2.count.csv")))
meta_data <- read.csv(paste0(cellbench_dir,"/RNAmix_celseq2.metadata.csv"))
all.equal(rownames(meta_data),colnames(rna_counts))
group <- apply(meta_data[,c('H2228_prop','H1975_prop','HCC827_prop')],1,function(x) paste0(x,collapse=" "))
H2228_to_H1975 <- H2228_to_HCC827 <- H1975_to_HCC827 <- rep(NA,length(group))

H2228_to_H1975[group=="1 0 0"] = 0
H2228_to_H1975[group=="0.68 0.16 0.16"] = 1
H2228_to_H1975[group=="0.33 0.33 0.33"] = 2
H2228_to_H1975[group=="0.16 0.68 0.16"] = 3
H2228_to_H1975[group=="0 1 0"] = 4

H2228_to_HCC827[group=="1 0 0"] = 0
H2228_to_HCC827[group=="0.68 0.16 0.16"] = 1
H2228_to_HCC827[group=="0.33 0.33 0.33"] = 2
H2228_to_HCC827[group=="0.16 0.16 0.68"] = 3
H2228_to_HCC827[group=="0 0 1"] = 4

H1975_to_HCC827[group=="0 1 0"] = 0
H1975_to_HCC827[group=="0.16 0.68 0.16"] = 1
H1975_to_HCC827[group=="0.33 0.33 0.33"] = 2
H1975_to_HCC827[group=="0.16 0.16 0.68"] = 3
H1975_to_HCC827[group=="0 0 1"] = 4
traj_orders <- data.frame(H2228_to_H1975,H2228_to_HCC827,H1975_to_HCC827)
#meta_data[which(apply(traj_orders,1,function(x) all(is.na(x)))),]

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"rnamixcs2_rnaseq_traj_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"rnamixcs2_rnaseq_traj_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"rnamixcs2_rnaseq_traj_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"rnamixcs2_rnaseq_traj_hvg500_expr.mtx"))
write.table(traj_orders,paste0(output_dir,"rnamixcs2_rnaseq_traj_labs.txt"),quote=F,row.names=F)

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

#cellbench sortseq RNA-mix
rna_counts <-as.matrix(read.csv(paste0(cellbench_dir,"/RNAmix_sortseq.count.csv")))
meta_data <- read.csv(paste0(cellbench_dir,"/RNAmix_sortseq.metadata.csv"))
all.equal(rownames(meta_data),colnames(rna_counts))
group <- apply(meta_data[,c('H2228_prop','H1975_prop','HCC827_prop')],1,function(x) paste0(x,collapse=" "))
H2228_to_H1975 <- H2228_to_HCC827 <- H1975_to_HCC827 <- rep(NA,length(group))

H2228_to_H1975[group=="1 0 0"] = 0
H2228_to_H1975[group=="0.68 0.16 0.16"] = 1
H2228_to_H1975[group=="0.33 0.33 0.33"] = 2
H2228_to_H1975[group=="0.16 0.68 0.16"] = 3
H2228_to_H1975[group=="0 1 0"] = 4

H2228_to_HCC827[group=="1 0 0"] = 0
H2228_to_HCC827[group=="0.68 0.16 0.16"] = 1
H2228_to_HCC827[group=="0.33 0.33 0.33"] = 2
H2228_to_HCC827[group=="0.16 0.16 0.68"] = 3
H2228_to_HCC827[group=="0 0 1"] = 4

H1975_to_HCC827[group=="0 1 0"] = 0
H1975_to_HCC827[group=="0.16 0.68 0.16"] = 1
H1975_to_HCC827[group=="0.33 0.33 0.33"] = 2
H1975_to_HCC827[group=="0.16 0.16 0.68"] = 3
H1975_to_HCC827[group=="0 0 1"] = 4
traj_orders <- data.frame(H2228_to_H1975,H2228_to_HCC827,H1975_to_HCC827)
#meta_data[which(apply(traj_orders,1,function(x) all(is.na(x)))),]

sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts,file=paste0(output_dir,"rnamixss_rnaseq_traj_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"rnamixss_rnaseq_traj_hvg5k_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"rnamixss_rnaseq_traj_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"rnamixss_rnaseq_traj_hvg500_expr.mtx"))
write.table(traj_orders,paste0(output_dir,"rnamixss_rnaseq_traj_labs.txt"),quote=F,row.names=F)

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))
