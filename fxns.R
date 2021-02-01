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
  diss <- diss[upper.tri(diss)]

  #s+, number of times intra-clsuter distances are strictly smaller than inter-cluster-distances
  sp <- sum(sapply(which(ind_vec), function(i) sum(diss_vec[which(!ind_vec)] > diss[i])))
  #s- number of times intra-clsuter distances are strictly greater than inter-cluster-distances   
  sp <- sum(sapply(which(ind_vec), function(i) sum(diss_vec[which(!ind_vec)] < diss[i])))
 
  tmp <- c(
    #calculation time
    time_tmp,
    #Baker-Hubert
    (sp - sm)/(sp + sm),
    #Gplus
    (2*sm)/(Nt*(Nt-1)),
    #Tau
    (sp - sm)/sqrt(Nb*Nw*((Nt*(Nt-1))/2))
  )

  return(c(
    "Calculation time (s)",
    "Baker-Hubert Gamma", = 
    "G-Plus",
    "Tau",
  ))

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


