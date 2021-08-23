#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
p_load(Matrix,RcppCNPy,reticulate,pbapply,cluster,purrr)
np <-import("numpy")

n_sim <- 30
n_cell <- 2830
set.seed(1354)

#calculate w for k=1 and k=2
calc_w <- function(dist_tmp,k){
  if(k==1){
    ut_tmp <- as.matrix(dist_tmp)
    ut_tmp <- ut_tmp[upper.tri(ut_tmp)]
    w_tmp <- sum(ut_tmp) / ( 2*attr(dist_tmp,"Size") )
  } else {
    ##induce PAM
    clust_tmp <- pam(x=dist_tmp,k=k,diss=T)$clustering
    nr_tmp <- as.numeric(table(clust_tmp))
    drs_tmp <- sapply(1:k, function(r) {
      idx_tmp <- which(clust_tmp==r)
      sbt_tmp <- as.matrix(dist_tmp)[idx_tmp,idx_tmp]
      sbt_tmp <- sbt_tmp[upper.tri(sbt_tmp)]
      dr_tmp <- sum(sbt_tmp)/( 2 * nr_tmp[r] )
    })
   w_tmp <- sum(drs_tmp)
  }
  return(w_tmp)
}

#function to simulate a unif(0,1) distance matrix
sim_randist <- function(n){
  nrv <- (n*(n-1))/2
  rval <- runif(n=nrv,min=0,max=1)
  tmp <- matrix(nrow=n,ncol=n,0)
  tmp[lower.tri(tmp,diag=F)] <- rval
  tmp <- t(tmp)
  tmp[lower.tri(tmp,diag=F)] <- rval
  tmp <- as.dist(tmp)
  return(tmp) 
}

#load in real distance, calculate W for that
k_tmp <- 5
dist_mat <- sim_randist(n_cell)
ws_true <- sapply(1:k_tmp, function(j) calc_w(dist_mat,j))

#loop over 1:nsim, calculate Ws for each
ws_bg <- sapply(1:n_sim, function(i) {
  dist_mat <- sim_randist(n_cell)
  ws_tmp <- sapply(1:k_tmp, function(j) calc_w(dist_mat,j))
  return(ws_tmp)
})

gap_stats <- sapply(1:k_tmp, function(k) mean(log(ws_bg[k,])) - log(ws_true[k]) )

outfile <- "/users/ndyjack/Dist_Proj/tables/gapstat/randdist.txt"
cat("\n writing output to ",outfile,"\n")
write(gap_stats,file=outfile)
