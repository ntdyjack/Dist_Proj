import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
import minocore
from scipy.io import mmread

args = sys.argv
#for x in args:
#  print(x)
#example args
#args = ['/users/ndyjack/Dist_Proj/scripts/calc_hash.py', '/users/ndyjack/Dist_Proj/tables/test_datasets/','/users/ndyjack/Dist_Proj/tables/hashes/','cellmix_rnaseq_traj_full','s2jsd']
mtx_dir = args[1] + args[3] + '_expr.mtx'
#ident_dir = args[1] + args[3] + '_labs.mtx'
output_dir = args[2] + args[3] + '.' + args[4] + '.csv'
seed = 1337

dat = mmread(mtx_dir) #read sparse matrix
dat = dat.toarray().transpose() #transpose and reformat
dat = dat / np.sum(dat, axis=1)[:,np.newaxis] #normalize
n_features = dat.shape[1] # get number of features

if args[4] == "cl1":
  hashf = minocore.ClippedL1LSHasher(dim=n_features, k=100, l=100, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif args[4] == "jsd":
  hashf = minocore.JSDLSHasher(dim=n_features, k=100, l=100, r=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif args[4] == "l1":
  hashf = minocore.L1LSHasher(dim=n_features, k=100, l=100, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif args[4] == "l2":
  hashf = minocore.L2LSHasher(dim=n_features, k=100, l=100, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif args[4] == "s2jsd":
  hashf = minocore.S2JSDLSHasher(dim=n_features, k=100, l=100, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
else:
  sys.exit("invalid hash function")

hashout_round = np.ceil(hashout).astype(np.uint16)
np.savetxt(output_dir, hashout_round, delimiter=",",fmt='%i')#store hash as csv to do computations
