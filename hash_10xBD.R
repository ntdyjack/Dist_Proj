args <- commandArgs(trailingOnly=T)
cat("\n","running Rscript","\n")
cat("args=",args)
library(Matrix)
library(DelayedArray)
library(HDF5Array)

#args= c("/users/ndyjack/Dist_Proj/tables/test_datasets/", "/fastscratch/myscratch/ndyjack/tables/", "/fastscratch/myscratch/ndyjack/hashes/", "400", "521494", "522800NULL")
args[6] <- sub("NULL","",args[6])

#example args
#data dir, mtx dir, hash dir, i, istart, i end 
# cells = 1306127
# genes = 24015

#setRealizationBackend("RleArray")
setRealizationBackend("HDF5Array")
dat <- HDF5Array(paste0(args[1],"1M_neurons_full_expr.h5"),name = "10xBD_Full")
dat <- as.matrix(realize(dat[,as.integer(args[5]):as.integer(args[6])]))
dat <- Matrix(dat,sparse=T)
writeMM(dat,paste0(args[2],"tmp_10xBD.",args[4],".mtx"))

#dat <- h5mread(paste0(args[1],"1M_neurons_full_expr.h5"),name="10xBD_Full")
#block_dim <- c(24015,as.integer(args[5]) - as.integer(args[4]))
#view_port <- ArrayViewport(as.integer(c(24015,1306127)),IRanges(c(1,as.integer(args[4])),width=block_dim))
#tickmarks=list(1:1000,as.integer(args[4]):as.integer(args[5])))
#ArbitraryArrayGrid(tickmarks=list(1:1000,as.integer(args[4]):as.integer(args[5])))
#dat <- read_block(x=paste0(args[1],"1M_neurons_full_expr.h5"),viewport = view_port)
#dat <- as.matrix(counts(dat[,as.integer(args[4]):as.integer(args[5])]))
#dat <- Matrix(dat,sparse=T)
#writeMM(dat,paste0(args[1],"tmp_10xBD.",args[3],".mtx"))
