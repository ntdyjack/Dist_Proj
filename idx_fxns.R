internal_indices <- function(mat,clust,distval) {
  require(Rfgc)
  #arguments:
    #mat - expression matrix, cells in rows, genes in columns. class = numeric matrix
    #clust - cluster identity vector for the cells (rows of mat). class = integer vector
    #distval - which Rfgc distance metric to calculate. class = integer

  #argument testing
  if(class(mat)!="matrix"){ stop("Error: mat must be a matrix") }
  if(!is.integer(clust)){ stop("Error: clust must be an integer vector") }
  if(length(clust)!=nrow(mat)){ stop("Error: dimension mismatch") }

  #unique cluster IDs
  N <- nrow(mat)
  cid_vec <- unique(clust)
  K <- length(cid_vec)
  cid_list <- lapply(cid_vec, function(i) which(clust==i) )
  cid_idx <- 1:K

  #index for the barrycenter vectors
  sum_idx <- (nrow(mat)+1):(nrow(mat)+K+1)
  #calculate the barycenters (IE, within-cluster means vectors)
  mu_overall <- colMeans(mat)
  muk_mat <- t(sapply(cid_list, function(i) colMeans(mat[i,])))
  mat <- rbind(mat,muk_mat,mu_overall)

  N <- nrow(mat)
  #Nt is the number of distinct pairs
  Nt <- length(clust)*(length(clust)-1)/2
  #Nw is the  total number of distinct pairs within each cluster
  Nw <- sum(sapply(cid_list, function(i) length(i)*(length(i)-1)/2))
  #Nb is the total numer of distinct pairs outside of each cluster 
  Nb <- Nt - Nw 

  #calculate and time distance 
  start_time <- Sys.time()
  diss <- Rfgc::dist_matrixdd(mat, distval-1)
  time_tmp <- Sys.time() - start_time 
  if(max(diss)>1){
    cat("scaling...\n")
    diss <- diss/max(diss)
  }

  #Sw is the sum of the within-cluster distances
  Sw <- sum(sapply(cid_list, function(i) sum(diss[i,i][upper.tri(diss[i,i])]) ))
  #Smin is the sum of the Nw smallest distances between all pairs of points
  Smin <- sum(sort(diss[-sum_idx,-sum_idx][upper.tri(diss[-sum_idx,-sum_idx])],decreasing=F)[1:Nw])
  #Smax is the sum of the Nw largest distances between all pairs of points
  Smax <- sum(sort(diss[-sum_idx,-sum_idx][upper.tri(diss[-sum_idx,-sum_idx])],decreasing=T)[1:Nw])
  #Sb is the sum of the between-cluster_distances
  Sb <- sum(sapply(cid_list, function(i) sum(diss[i,-c(i,sum_idx)][upper.tri(diss[i,-c(i,sum_idx)])])))
  #BGSS is the weighted sum of the squared distances between each barrycenter the overall barrycenter
  BGSS <- sum(sapply(cid_idx, function(i) length(cid_list[[i]])*diss[sum_idx[i],sum_idx[K+1]]^2 )) 
  #WGSS is the within cluster dispersion (within-cluster sum of squared distance)
  WGSS <- sum(sapply(cid_idx, function(i) sum((diss[cid_idx[[i]],sum_idx[i]]^2)) ))
  #dk is the  mean distance of the points belonging to cluster to thir barrycenter
  dk <- sapply(cid_idx, function(i) mean((diss[cid_idx[[i]],sum_idx[i]]^2)) ) 
  #dmin smallest distance between 2 differnet clusters
  dmin <- min(sapply(cid_list, function(i) min(diss[i,-c(i,sum_idx)])))
  #dmax largest difference within a cluster
  dmax <- max(sapply(cid_list, function(i) max(diss[i,i])))
  #Mb: the smallest distance value between 2 points not in the same cluster
  #Mb <- min(sapply(cid_list, function(i) min(diss[i,-c(i,sum_idx)])))
  #s+, number of times intra-clsuter distances are strictly smaller than inter-cluster-distances
  sp <- sum(sapply(cid_list, function(i) { 
    sum(sapply(diss[i,i][upper.tri(diss[i,i])], function(x) {
      sum(x < diss[i,-c(i,sum_idx)][upper.tri(diss[i,-c(i,sum_idx)])])
    }))
  })) 
  #sum(sapply(cid_list, function(i)  sum(diss[i,i][upper.tri(diss[i,i])] < Mb)))
  #s- number of times intra-clsuter distances are strictly greater than inter-cluster-distances
  sm <- sum(sapply(cid_list, function(i) {
    sum(sapply(diss[i,i][upper.tri(diss[i,i])], function(x) {
      sum(x > diss[i,-c(i,sum_idx)][upper.tri(diss[i,-c(i,sum_idx)])])
    }))
  })) 
  #sum(sapply(diss[i,i][upper.tri(diss[i,i])], function(x) sum(x > diss[i,-c(i,sum_idx)][upper.tri(diss[i,-c(i,sum_idx)])])))
  #sum(sapply(cid_list, function(i)  sum(diss[i,i][upper.tri(diss[i,i])] > Mb)))
  #Db: largest distance between two barrycenters
  Db <- max(diss[sum_idx[-(K+1)],sum_idx[cid_idx]])
  #Dm: smallest distance between two barrcenters
  Dm <- min(diss[sum_idx[cid_idx],sum_idx[cid_idx]][upper.tri(diss[sum_idx[cid_idx],sum_idx[cid_idx]])])
  #Ew sum of the distances of points to their barycenter
  Ew <- sum(sapply(cid_idx, function(i) sum(diss[cid_list[[i]],sum_idx[i]]))) 
  #Et sum of the distances of points to the overall barycenter 
  Et <- sum(diss[-sum_idx,sum_idx[K+1]]) 
  #sum(sapply(cid_idx, function(i) sum(diss[cid_list[[i]],sum_idx[K+1]]))) 
  #diss[cid_list[[i]],sum_idx[K+1]][upper.tri(diss[cid_list[[i]],sum_idx[K+1]])] )))
  #Rm: mean of quotient of the distance from one point to its barycenter / smallest distance to another barycenter
  #Rm <-
  #Jk <-  sapply(rm, function(i) max(c(0,1)-i))
 
  tmp <- c(
    #calculation time
    time_tmp,
    #Ball-Hall: mean of the mean of squared distance from each point to its barycenter
    mean(sapply(cid_idx, function(i) mean((diss[cid_list[[i]],sum_idx[i]]^2)) )), 
    #C-Index: 
    (Sw - Smin) / (Smax - Smin),
    #Calinski-Harabasz:
    ((N-K)*(BGSS)) / ((WGSS)*(K-1)),
    #Davies-Bouldin
    mean(sapply(cid_idx, function(i) max((dk[-i]+dk[i])/(diss[sum_idx[i],sum_idx[-c(i,K+1)]])) )), 
    #Dunn
    dmin/dmax, 
    #Baker-Hubert
    (sp - sm)/(sp + sm),
    #Generalized Dunn
    #### To-do     
    #Gplus
    (2*sm)/(Nt*(Nt-1)),
    #Log SS Ratio
    log(BGSS/WGSS),
    #McClain-Rao:
    (Nb * Sw) / (Nw * Sb),
    #PBM
    ((Db*Et)/(K*Ew))^2,
    #Point-Biserial
    ((Sw/Nw) - (Sb/Nb)) * sqrt(Nw*Nb)/Nt,
    #Ray-Turi
    (WGSS/N) * (1/(Dm^2))#, 
    #Silhouette
    ####mean(sapply( ))
    #Tau
    ####(sp - sm)/sqrt(Nb*Nw*((Nt*(Nt-1)/2)))
    #Wammert-Gancarski
    ###mean(sapply(cid_idx, function(i) length(cid_list[[i]])*Jk[i]))#, 
    #Xie-Beni
    ####(1/N)*() 
  )

names(tmp) <- c(
    "Calculation time (s)",
    "Ball-Hall",
    "C-Index",
    "Calinski-Harabasz",
    "Davies-Bouldin",
    "Dunn",
    "Baker-Hubert",
#    "Generalized Dunn",
    "G-Plus",
    "Log SS Ratio",
    "McClain-Rao",
    "PBM",
    "Point Biserial",
    "Ray-Turi"#,
#    "Silhouette",
#    "Tau",
#    "Wammert-Gancarski",
#    "Xie-Beni"
  )

  return(tmp)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

draw_legend <- function(points=T){
  par(new = "TRUE",plt = c(0.03,0.97,0.03,0.23),las = 1, cex.axis = 1)
  plot(x=1,y=1,type='n',main="",xlab="",ylab="",xaxt="n", yaxt="n",bty='n')
  if(points){
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=4,
      x.intersp=0.6,y.intersp=1,text.width=0.2,cex=1.0,
      pch=shps,pt.bg="#00000050",pt.cex=1.0,bty='n')
  }else{
    legend('left',legend=dists_vec[idx_test],col=cols,ncol=4,
      x.intersp=0.5,y.intersp=1,text.width=0.2,cex=1.0,
      lty=ltys,lwd=1.5,bty='n')
  }
}


