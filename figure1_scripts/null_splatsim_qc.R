#ml old_conda_R/3.5.x
options("width"=180)
library(pacman)
p_load(DropletUtils,heatmap3,viridis,readxl,Matrix,scran,splatter,mbkmeans,scater,scry,purrr)

output_dir <- "/users/ndyjack/Dist_Proj/tables/test_datasets/"

#############################
######simulate 1000 splatter cells with no cluster structure
#############################

cex_min <- 0.5
cex_max <- 2.5
map_cex <- function(x,y1,y2) {
  cex_min + ((cex_max - cex_min) / (y2 - y1) * (x - y1 ) )
}

#splatter simulation paramters
seed <- 1472
#number of genes
nG <- 20000
#n, total number of cells
nC <- 1000
#lS library size
lS <- 9.5

#initial parameters
params <- newSplatParams()
params <- setParam(params, "nGenes" , nG)
params <- setParam(params, "seed", seed)
params <- setParam(params, "lib.loc", lS)
params <- setParam(params, "batchCells" , nC)

dat <- splatSimulate(params ,verbose = FALSE)
rna_counts <- as.matrix(dat@assays@data$counts)
drop <- which(rowSums(rna_counts)==0)
if(length(drop)>0){ rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),] }
outlier <- rowData(dat)[,'OutlierFactor']!=1
libsizes <- colSums(rna_counts)
lim <- c(min(libsizes),max(libsizes))
ptsizes <- sapply(libsizes, function(z) map_cex(z,lim[1],lim[2]))
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1),counts=rna_counts))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
sce <- computeSumFactors(sce, min.mean=0.1)
sce <- logNormCounts(sce)
sce <- scater::runPCA(sce, ncomponents = 50,subset_row=hvg_1k,scale = TRUE, name='hvg')

pdf("/users/ndyjack/Dist_Proj/images/nullsim_pcaplot.pdf",width=6,height=6)
  plot(x= reducedDims(sce)[[1]][,1],y= reducedDims(sce)[[1]][,2],
    main="",xlab="",ylab="",xaxt="n", yaxt="n",
    col='black',pch=ifelse(outlier,13,16),cex=ptsizes)
dev.off()

#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"splatsim_1cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"splatsim_1cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"splatsim_1cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"splatsim_1cl_hvg5k_expr.mtx"))

#simulate background uniform datasets for the ARI tests
test_prefix <- "/users/ndyjack/Dist_Proj/tables/null_datasets/splat_bgsim_1cl"
feat_ranges <- apply(rna_counts,1,function(x) c(min(x),max(x)))
n_sim <- 30
n <- ncol(rna_counts)
m <- nrow(rna_counts)
set.seed(1234)
for(j in 1:30){
  name_tmp <- paste0(test_prefix,"_",j)
  dat_tmp <- t(sapply(1:m, function(i)  rdunif(n=n,a=feat_ranges[1,i],b=feat_ranges[2,i])))
  rownames(dat_tmp) <- rownames(rna_counts)
  dat_tmp <- Matrix(dat_tmp,sparse=T)
  writeMM(obj = dat_tmp, file=paste0(name_tmp,"_full_expr.mtx"))
  writeMM(obj = dat_tmp[hvg_500,], file=paste0(name_tmp,"_hvg500_expr.mtx"))
  writeMM(obj = dat_tmp[hvg_1k,], file=paste0(name_tmp,"_hvg1k_expr.mtx"))
  writeMM(obj = dat_tmp[hvg_5k,], file=paste0(name_tmp,"_hvg5k_expr.mtx"))
}


