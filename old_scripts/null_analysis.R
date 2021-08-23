#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
p_load(Matrix,RcppCNPy,reticulate,pbapply,cluster,purrr)
#np <-import("numpy")

#workflow....

#Wk = sum_r D_r / (2*n_r)
#Gap(k) = (1/B) * sum_b { log(W_kb) - log(W_k) }


#llr_fxns
llr <- function(counts1,counts2){
  nz1 <- counts1>0
  nz2 <- counts2>0
  nze <- counts1 | counts2
  prob1 <- counts1 / sum(counts1[nz1])
  prob2 <- counts2 / sum(counts2[nz2])
  llha1 <- counts1[nz1] * log(prob1[nz1])
  llha2 <- counts2[nz2] * log(prob2[nz2])
  llha <- sum(llha1) + sum(llha2)
  prob_avg <- (prob1 + prob2) / 2
  lg_pavg <- log(prob_avg[nze])
  llhn1 <- counts1[nze] * lg_pavg
  llhn2 <- counts2[nze] * lg_pavg
  llhn <- sum(llhn1) + sum(llhn2)
  return(llha - llhn)
}

compd <- function(dat){
  n <- ncol(dat)
  dist_mat <- matrix(0,nrow=n,ncol=n)
  dist_vec <- unlist(lapply(1:(n-1), function(i) {
    sapply((i+1):(n), function(j) {
      d <- llr(dat[,i],dat[,j])
      return(d)
    })
  }))
  dist_mat[lower.tri(dist_mat)] <- dist_vec
  dist_mat <- t(dist_mat)
  dist_mat[lower.tri(dist_mat)] <- dist_vec
  return(as.dist(dist_mat))
}


#function that does the following:
calc_w <- function(x,flag){
  ##calculate distance
  if(flag) {
    dist_tmp <- dist(x)
  }else{
    dist_tmp <- compd(t(x))
  }
  ut_tmp <- as.matrix(dist_tmp)
  ut_tmp <- ut_tmp[upper.tri(ut_tmp)]
  ##calculate W1 (just the sum of all distnaces)
  w1_tmp <- sum(ut_tmp) / ( 2*nrow(x) )
  ##induce PAM
  clust_tmp <- pam(x=dist_tmp,k=2,diss=T)$clustering
  nr_tmp <- as.numeric(table(clust_tmp))
  ##calculate W2 (sum of within-cluster differences)
  drs_tmp <- sapply(1:2, function(r) {
    idx_tmp <- which(clust_tmp==r)
    sbt_tmp <- as.matrix(dist_tmp)[idx_tmp,idx_tmp]
    sbt_tmp <- sbt_tmp[upper.tri(sbt_tmp)]
    dr_tmp <- sum(sbt_tmp)/( 2 * nr_tmp[r] )
  })
  w2_tmp <- sum(drs_tmp)
  ##return w1 and w2
  return(c(w1_tmp,w2_tmp)) 
}

#read in null dataset (full)
dat <- readMM('/users/ndyjack/Dist_Proj/tables/test_datasets/only293t_rnaseq_1cl_hvg500_expr.mtx') 
#readMM('/users/ndyjack/Dist_Proj/tables/test_datasets/splat_sim_1cl_hvg500_expr.mtx')
#number cells to do set on
n <- 100
idx <- 1:n
dat <- as.matrix(dat[,idx])
#dat <- dat[-which(rowSums(dat)==0),]
#get feature ranges for each
feat_ranges <- apply(dat,1,function(x) c(min(x),max(x)))
dat <- t(dat)
#calculate initial w1 and w2
true_ws_l2 <- calc_w(dat,T)
true_ws_lr <- calc_w(dat,F)
#set seed, loop over 1:nsim (30?)
set.seed(123)
n_sim <- 30
bg_ws <- pblapply(1:n_sim, function(b) {
  dat_tmp <- sapply(1:ncol(dat), function(i)  rdunif(n=n,a=feat_ranges[1,i],b=feat_ranges[2,i]))
  ws_tmp_l2 <- calc_w(dat_tmp,T)
  ws_tmp_lr <- calc_w(dat_tmp,F)
  list(l2=ws_tmp_l2,lr=ws_tmp_lr)
})

bg_ws_l2 <- sapply(bg_ws, function(x) x$l2)
bg_ws_lr <- sapply(bg_ws, function(x) x$lr)

gap_stats_l2 <- sapply(1:2, function(k) mean(log(bg_ws_l2[k,])) - log(true_ws_l2[k]) )
gap_stats_lr <- sapply(1:2, function(k) mean(log(bg_ws_lr[k,])) - log(true_ws_lr[k]) )

gap_stats_lr[1]/gap_stats_lr[2]
gap_stats_l2[1]/gap_stats_l2[2]

#simulate a discrete uniform (n=100) with the same row ranges as original dataset
#calculate w1 and w2 with the calc_w function

#d_l2 <- np$load('/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_full.1.npy')
#d_l2 <- d_l2[idx,idx]
#n <- n
#d_l2 <- d_l2[upper.tri(d_l2)]
#d_bh <- np$load('/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_full.7.npy')
#d_bh <- d_bh[idx,idx]
#d_bh <- d_bh[upper.tri(d_bh)]


#n_sim <- 1000

#set.seed(12)
#for each simulation, make a random
#sim_gp <- pbsapply(1:n_sim, function(i) {
  #simulat random labels rounding from unif(0,1)
#  cluster_lab <- runif(n=n,min=0,max=1)
#  cluster_lab <- round(cluster_lab)
#  ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
#  tmp <- ind_vec[upper.tri(ind_vec,diag=F)]
#  ind_intra <- which(tmp)
#  ind_inter <- which(!tmp)

  #estimate g+
#  n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
#  Nt_samp <- (n_samp + n_samp)*(n_samp + n_samp - 1)/2

#  gp_cur <- sapply(list(d_l2,d_bh), function(dist_vec) {
#    dist_inter <- sort(dist_vec[ind_inter])
#    dist_intra  <- sort(dist_vec[ind_intra])
#    orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
#    orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
#    gp <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nt_samp
#    return(gp)
#  })
#
#  return(gp_cur)
#})

#args <- c("/users/ndyjack/Dist_Proj/rdata/","/users/ndyjack/Dist_Proj/tables/","jurkat293t_rnaseq_2cl_hvg500","0","gplus_results_July2020.txt")
#cell cluster identity vector
#get rid of hvg string
#cluster_lab <- as.integer(readLines(paste0(args[2],sub("_[h|f].*","",args[3]),"_labs.txt")))


#read in cell distance vector
#dist_mat
#dist_vec <- readRDS(paste0(args[1],args[3],".",args[4],".rda"))$'dist'

#load in L2


#load in Bhatacharrya

#g+ parameters
#Nt <- length(cluster_lab)*(length(cluster_lab)-1)/2
#Nz <- (Nt*(Nt-1))/2
#cell pair x cluster identity vector
#ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
#ind_vec <- apply(ind_vec,1,function(x) {
#  tmp <- rep(0,length(x))
#  tmp[which(x)] <- cluster_lab[which(x)]
#  return(tmp)
#})
#ind_vec <- ind_vec[upper.tri(ind_vec)]
#ind_intra <- which(ind_vec>0)
#ind_inter <- which(ind_vec==0)

#cat("\nestimating g+...\n")
#n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
#Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
#dist_inter <- sort(dist_vec[ind_inter])
#dist_intra  <- sort(dist_vec[ind_intra])
#rm(dist_vec) #memory cleanup
#orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
#orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
#gplus_samp <-  sum(sapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nz_samp
#cat("\nfinished g+...\n")

#output <- paste0(args[3],".",args[4]," ",gplus_samp)
#outfile <- paste0(sub(x=args[2],pattern="test_datasets","gplus"),args[3],".",args[4],".txt")
#cat("\n writing output to ",outfile,"\n")
#write(output,file=outfile)
