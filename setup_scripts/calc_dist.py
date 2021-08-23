import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
import minicore as mc
from scipy.io import mmread

args = sys.argv
for x in args:
  print(x)
#example args
#args = ['/users/ndyjack/Dist_Proj/scripts/calc_hash.py', '/users/ndyjack/Dist_Proj/tables/test_datasets/','/users/ndyjack/Dist_Proj/pickles/','jurkat293t_rnaseq_2cl_hvg500','0']
#mtx_dir = '/users/ndyjack/Dist_Proj/tables/null_datasets/splat_bgsim_1cl_10_hvg500_expr.mtx'
mtx_dir = args[1] + args[3]# + '_expr.mtx'
output_dir = args[2] + args[3] + '.' + args[4] + '.npy'

dat = mmread(mtx_dir).transpose() #read sparse matrix
dat = dat.tocsr() #convert to row sparse matrix

counts = np.asarray(dat.data,dtype=np.uint16) #non-zero values
idx = np.asarray(dat.indices,dtype=np.uint16) #indices
ip = np.asarray(dat.indptr,dtype=np.uint32) #index pointers
shape = np.asarray(dat.shape,dtype=np.uint32) #shape
nonzer = dat.getnnz()

dat = mc.csr_tuple(data=counts,indices=idx,indptr=ip,shape=shape,nnz=nonzer) #set up csr tuple
dat = mc.CSparseMatrix(dat) #set up CSparseMatrix

#check if we need prior
m = int(args[4])

if m in [5,13,14,23]:
  p = 1
else:
  p = 0


dist = mc.cmp(matrix = dat, data = dat,msr = m, prior = p)  #compute distance
np.save(output_dir,dist) #save as pickle
