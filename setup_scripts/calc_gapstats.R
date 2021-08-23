#load accessory functions, set parameters
args <- commandArgs(trailingOnly=T)
cat("\n",args,"\n")
if(!require(pacman)){ install.packages("pacman") }
p_load(Matrix,RcppCNPy,reticulate,pbapply,cluster,purrr)
np <-import("numpy")

n_sim <- 30
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

#args <-c("/users/ndyjack/Dist_Proj/pickles/","/fastscratch/myscratch/ndyjack/pickles/","/users/ndyjack/Dist_Proj/tables/gapstat/","hvg1k","10")

#name string?
#dt_str <- paste0(args[1],"only293t_rnaseq_1cl_",args[4],".",args[5],".npy")
dt_str <- paste0(args[2],"splatsim_1cl_",args[4],".",args[5],".npy")
#lab_dir <- sub("pickles","tables/test_datasets",dt_str)
#lab_dir <- sub("npy","labs.txt",lab_dir)
k_tmp <- 5

#load in real distance, calculate W for that
#dist_mat <- np$load("/users/ndyjack/Dist_Proj/pickles/splat_sim_1cl_hvg5k.13.npy")
#dist_mat <- np$load("/fastscratch/myscratch/ndyjack/pickles/splat_bgsim_1cl_1_full.13.npy")
dist_mat <- np$load(dt_str)
dist_mat <- as.dist(dist_mat)
ws_true <- sapply(1:k_tmp, function(j) calc_w(dist_mat,j))

#loop over 1:nsim, calculate Ws for each
ws_bg <- sapply(1:n_sim, function(i) {
  #bg293t_sim_1cl_14_hvg5k.7.npy
  #fname <- paste0(args[2],"bg293t_sim_1cl_",i,"_",args[4],".",args[5])
  fname <- paste0(args[2],'splat_bgsim_1cl_',i,"_",args[4],".",args[5])
  bg_dir <- paste0(fname,".npy")
  #lab_dir <- paste0(fname,".txt")
  #lab_dir <- sub("pickles","labels")
  dist_mat <- np$load(bg_dir)
  dist_mat <- as.dist(dist_mat)
  ws_tmp <- sapply(1:k_tmp, function(j) calc_w(dist_mat,j))
  return(ws_tmp)
})

gap_stats <- sapply(1:k_tmp, function(k) mean(log(ws_bg[k,])) - log(ws_true[k]) )
#gsr <- gap_stats[1]/gap_stats[2]

#output <- paste0(args[3],".",args[4]," ",km_ari)
#outfile <- paste0(args[3],'',args[4],".",args[5],".txt")
outfile <- paste0(args[3],'splat.',args[4],".",args[5],".txt")
cat("\n writing output to ",outfile,"\n")
write(gap_stats,file=outfile)

#calculate GS for both
#calcualte GSR
#write output

#read in null dataset (full)
#dat <- readMM('/users/ndyjack/Dist_Proj/tables/test_datasets/only293t_rnaseq_1cl_hvg500_expr.mtx') 
#readMM('/users/ndyjack/Dist_Proj/tables/test_datasets/splat_sim_1cl_hvg500_expr.mtx')
#number cells to do set on
#n_sim <- 30
#bg_ws <- pblapply(1:n_sim, function(b) {
#  dat_tmp <- sapply(1:ncol(dat), function(i)  rdunif(n=n,a=feat_ranges[1,i],b=feat_ranges[2,i]))
#  ws_tmp_l2 <- calc_w(dat_tmp,T)
#  ws_tmp_lr <- calc_w(dat_tmp,F)
#  list(l2=ws_tmp_l2,lr=ws_tmp_lr)
#})

#bg_ws_l2 <- sapply(bg_ws, function(x) x$l2)
#bg_ws_lr <- sapply(bg_ws, function(x) x$lr)
#gap_stats_l2 <- sapply(1:2, function(k) mean(log(bg_ws_l2[k,])) - log(true_ws_l2[k]) )
#gap_stats_lr <- sapply(1:2, function(k) mean(log(bg_ws_lr[k,])) - log(true_ws_lr[k]) )
#gap_stats_lr[1]/gap_stats_lr[2]
#gap_stats_l2[1]/gap_stats_l2[2]


#output <- paste0(args[3],".",args[4]," ",km_ari)
#outfile <- paste0(sub(x=args[2],pattern="test_datasets","gsr"),args[3],".",args[4],".kmeds.txt")
#cat("\n writing output to ",outfile,"\n")
#write(output,file=outfile)
