#script to submit a job for each test data-set
import os
import sys
import time
import numpy as np

#idx_test = ['cl1','jsd','l1','l2','s2jsd']
sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_hashBD.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/hash_10xBD.R '
py_script = '/users/ndyjack/Dist_Proj/scripts/hash_10xBD.py '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
data_dir = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '
mtx_dir = "/fastscratch/myscratch/ndyjack/tables/ "
hash_dir = "/fastscratch/myscratch/ndyjack/hashes/ "

template_arg = sh_script + r_script + py_script + data_dir + mtx_dir + hash_dir

n_cells = 1306127
n_jobs = 1000
l = int(np.ceil(n_cells/n_jobs))


for i in range(1,n_jobs+1):
  idx_use = [(i-1)*l+1,i*l]
  if i==n_jobs:
    idx_use[1] = n_cells

  args = template_arg + str(i) + ' ' + str(idx_use[0]) + ' ' + str(idx_use[1])
  name = 'hash10xBD' + '.' + str(i)
  err = ' -e ' + out_dir + name  + '.err '
  out = '-o ' + out_dir + name + '.out '
  cmd = 'qsub -N ' + name + err + out + args
  #print(cmd)
  os.system(cmd)
  time.sleep(1)
