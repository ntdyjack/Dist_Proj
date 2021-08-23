#ml old_conda_R/3.5.x
options("width"=180)
library(pacman)
p_load(DropletUtils,heatmap3,viridis,readxl,Matrix,scran,splatter,mbkmeans,scater,scry,purrr,scDesign2)

output_dir <- "/users/ndyjack/Dist_Proj/tables/test_datasets/"

#############################
#######10x 293t only data
#############################

#load in data from jurkat 293t 50/50
sce <- read10xCounts("/users/ndyjack/Dist_Proj/tables/293t_only/filtered_matrices_mex/hg19")
#make identity vector for jurkat vs 239t
rna_counts <- as.matrix(counts(sce))
rna_colsums <- colSums(rna_counts[grep("ENSG",rownames(rna_counts)),],)
mito_genes <- rownames(rna_counts)[grep("^MT-",sce@rowRanges@elementMetadata@listData$Symbol)]
pct_mito <- colSums(rna_counts[mito_genes,]) / rna_colsums
#remove zero-expression cite-seq cell
qtls_dropm <- quantile(pct_mito,c(0.005,0.995))
qtls_dropr <- quantile(rna_colsums,c(0.005,.995))
cells_drop <- which(
  rna_colsums > qtls_dropr[2] |  rna_colsums < qtls_dropr[1] |
  pct_mito > qtls_dropm[2] |  pct_mito < qtls_dropm[1]
)


pdf("/users/ndyjack/Dist_Proj/images/only293t_qc.pdf",width=12,height=8)
hist(breaks=100,rna_colsums)
abline(v=qtls_dropr,col="red",lty=2)
hist(breaks=100,pct_mito)
abline(v=qtls_dropm,col="red",lty=2)
dev.off()

#subset to cells that pass QC, normalize
rna_counts <- rna_counts[,-cells_drop]

#compute variant genes
rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),]
sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1),counts=rna_counts))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)

#normalize and do pca plot
sce <- computeSumFactors(sce, min.mean=0.1)
sce <- logNormCounts(sce)
sce <- scater::runPCA(sce, ncomponents = 50,ntop = 1000,scale = TRUE,BSPARAM = BiocSingular::RandomParam())
sce <- nullResiduals(sce, assay="counts", type="deviance")
sce <- scater::runPCA(sce, ncomponents = 50,ntop = 1000,exprs_values = "binomial_deviance_residuals",scale = TRUE, name = "GLM-PCA",BSPARAM = BiocSingular::RandomParam())

#colData(sce,'UMI') <- colSums(rna_counts)
sce@colData@listData$nUMI <- colSums(rna_counts) 

pdf("/users/ndyjack/Dist_Proj/images/only293t_pcaplots.pdf")
#plotPCA(sce)
plotReducedDim(sce, dimred = "PCA",colour_by='sizeFactor')
plotReducedDim(sce, dimred = "PCA",colour_by='nUMI')
plotReducedDim(sce, dimred = "GLM-PCA",colour_by='sizeFactor')
plotReducedDim(sce, dimred = "GLM-PCA",colour_by='nUMI')
dev.off()

#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"only293t_rnaseq_1cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"only293t_rnaseq_1cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"only293t_rnaseq_1cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"only293t_rnaseq_1cl_hvg5k_expr.mtx"))

#cat(dim(rna_counts),"\n")
#cat(length(clust),"\n")
#summary(colSums(rna_counts))
#summary(rowSums(rna_counts))

#simulate background uniform datasets for the ARI tests
test_prefix <- "/users/ndyjack/Dist_Proj/tables/null_datasets/bg293t_sim_1cl"
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


#############################
#######splatter-simulated data
#############################











#############################
#######splatter-simulated data
#############################

#splatter simulation paramters
seed <- 1472
#total number of cells
nC <- 1000
#probability of being in each group
#gP <- c(1.0)
#probability of difexp genes
#dP <- 0.0
#number of genes
nG <- 20000
#dropout shape (lower values mean more dropout if dO is positive)
dS <- -1
#library size factor
lS <- 10
#dropout midpoint
dO <- 1
#logfoldchange factor
lF <- 0

params <- newSplatParams()
params <- setParam(params, "nGenes" , nG)
params <- setParam(params, "seed", seed)
params <- setParam(params, "batchCells" , nC)
#params <- setParam(params, "group.prob", gP)
#params <- setParam(params, "de.prob", dP)
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.shape",dS)
params <- setParam(params, "lib.loc", lS)
params <- setParam(params, "dropout.mid", dO)
params <- setParam(params, "de.facLoc", lF)

#simulate data
dat <-  splatSimulate(params,method="single",verbose=F)
rna_counts <- as.matrix(dat@assays@data$counts)

pct_dropout <- formatC( 100*(sum(dat@assays@data$Dropout) / (nC*nG)),format='f',digits=1)
pct_zero <- formatC( 100*(sum(dat@assays@data$counts==0) / (nC*nG)),format='f',digits=1)

rna_colsums <- colSums(rna_counts)
qtls_dropr <- quantile(rna_colsums,c(0.005,.995))
cells_drop <- which( rna_colsums > qtls_dropr[2] |  rna_colsums < qtls_dropr[1] )

pdf("/users/ndyjack/Dist_Proj/images/splatter_qc.pdf",width=12,height=8)
hist(breaks=100,rna_colsums)
abline(v=qtls_dropr,col="red",lty=2)
dev.off()

rna_counts <- rna_counts[,-cells_drop]
rna_counts <- rna_counts[-which(rowSums(rna_counts)==0),]


sce <- SingleCellExperiment(assays=list(logcounts=log2(rna_counts+1)))
fit <- modelGeneVar(sce)
hvg_500 <- getTopHVGs(fit, n=500)
hvg_1k <- getTopHVGs(fit, n=1000)
hvg_5k <- getTopHVGs(fit, n=5000)
#write outputs to sparse matrix
rna_counts <- Matrix(rna_counts,sparse=T)
writeMM(obj = rna_counts, file=paste0(output_dir,"splat_sim_1cl_full_expr.mtx"))
writeMM(obj = rna_counts[hvg_500,], file=paste0(output_dir,"splat_sim_1cl_hvg500_expr.mtx"))
writeMM(obj = rna_counts[hvg_1k,], file=paste0(output_dir,"splat_sim_1cl_hvg1k_expr.mtx"))
writeMM(obj = rna_counts[hvg_5k,], file=paste0(output_dir,"splat_sim_1cl_hvg5k_expr.mtx"))

cat(dim(rna_counts),"\n")
cat(length(clust),"\n")
summary(colSums(rna_counts))
summary(rowSums(rna_counts))

