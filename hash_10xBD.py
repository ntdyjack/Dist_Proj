print("\n running hash")
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
import minocore
from scipy.io import mmread

args = sys.argv
print("args=",args)
#args = ["dontuse", "/users/ndyjack/Dist_Proj/tables/test_datasets/", "/users/ndyjack/Dist_Proj/tables/hashes/", "1000", "1305694", "1306127"]
#mtx_dir = '/fastscratch/myscratch/ndyjack/tables/tmp_10xBD.4.mtx'
mtx_dir = args[2] + 'tmp_10xBD.' + args[4] + '.mtx'
output_dir = args[3] + 'tmp_10xBD.' + args[4] + '.csv'
seed = 1337

dat = mmread(mtx_dir) #read sparse matrix
dat = dat.toarray().transpose() #transpose and reformat
dat = dat / np.sum(dat, axis=1)[:,np.newaxis] #normalize
n_features = dat.shape[1] # get number of features

#hasher=args[5]
hasher='jsd'

if hasher == "cl1":
  hashf = minocore.ClippedL1LSHasher(dim=n_features, k=10, l=10, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif hasher == "jsd":
  hashf = minocore.JSDLSHasher(dim=n_features, k=10, l=10, r=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif hasher == "l1":
  hashf = minocore.L1LSHasher(dim=n_features, k=10, l=10, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif hasher == "l2":
  hashf = minocore.L2LSHasher(dim=n_features, k=10, l=10, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
elif hasher == "s2jsd":
  hashf = minocore.S2JSDLSHasher(dim=n_features, k=10, l=10, w=0.01,seed=seed) #set up hash function
  hashout = hashf.project(dat) #project data
else:
  sys.exit("invalid hash function")

hashout_round = np.floor(hashout).astype(np.uint16)

tmp = np.sum(hashout_round,1)
p = [i for i, e in enumerate(tmp) if e == 0]

#[i for i, e in tmp if e == 0]
# [item for item in tmp if item == 0]
#np.sum(hashout_round,1) == 0


np.savetxt(output_dir, hashout_round, delimiter=",",fmt='%i')#store hash as csv to do computations

#for x in np.sum(hashout_round,1):
#  print(x)

#delete the temporary mtx file
