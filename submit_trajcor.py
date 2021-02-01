#script to submit a job for each test data-set
import os
import sys
import time


files = [
  "cellmix_rnaseq_traj",
  "rnamixcs2_rnaseq_traj",
  "rnamixss_rnaseq_traj"
]

gsets = [
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
]

fnames = [x+y for x in files for y in gsets]

#idx_test = [ [x for x in range(0,22)], [x for x in range(28,30)],["std_pca", "trn_pca", "glm_pca"]]
#idx_test = [item for sublist in idx_test for item in sublist]
idx_test = [5,6,12,13]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/calc_trajcor.R '
data_dir  = '/users/ndyjack/Dist_Proj/rdata/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '

for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i) + ' '
    name_tmp = x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.gplus.err '
    out = '-o ' + out_dir + name_tmp + '.gplus.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    #print(args)
    os.system(cmd)
    time.sleep(1)
