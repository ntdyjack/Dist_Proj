#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
p_load(Matrix,RcppCNPy,reticulate,pbapply)
np <-import("numpy")

d_l2 <- np$load('/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_hvg500.1.npy')
d_l2 <- d_l2[1:500,1:500]
n <- ncol(d_l2)
d_l2 <- d_l2[upper.tri(d_l2)]


#d_bh <- np$load('/users/ndyjack/Dist_Proj/pickles/only293t_rnaseq_1cl_hvg500.7.npy')
#d_bh <- d_bh[upper.tri(d_bh)]


#test if our estimator is actually right?
Nt <- n*(n-1)/2
Nz <- (Nt*(Nt-1))/2

#set.seed(12)
cluster_lab <- runif(n=n,min=0,max=1)
cluster_lab <- round(cluster_lab)
ind_vec <- sapply(cluster_lab, function(x) sapply(cluster_lab, function(y) x==y))
tmp <- ind_vec[upper.tri(ind_vec,diag=F)]
ind_intra <- which(tmp)
ind_inter <- which(!tmp)
gp_true <- sum(pbsapply(d_l2[ind_intra], function(x) sum(x > d_l2[ind_inter]))) / Nz

#what if we used a denominator of just the number of inter-intra comparisons...?
sm_true <- sum(pbsapply(d_l2[ind_intra], function(x) sum(x > d_l2[ind_inter])))
denom_test <- as.numeric(length(ind_intra)) * as.numeric(length(ind_inter))



n_samp <- round(0.005*length(ind_vec)) #sample 1% of all distances (0.5% inter, 0.5% intra)
Nt_samp <- (n_samp + n_samp)*(n_samp + n_samp - 1)/2
#Nz_samp <- (Nt_samp)*(Nt_samp-1)/2
dist_inter <- sort(d_l2[ind_inter])
dist_intra  <- sort(d_l2[ind_intra])
orders_inter <- round(seq(1,length(ind_inter),length.out=n_samp))
orders_intra <- round(seq(1,length(ind_intra),length.out=n_samp))
gp_samp <-  sum(pbsapply(dist_intra[orders_intra], function(x) sum(x > dist_inter[orders_inter]))) / Nt_samp

#test this for the estimator?

