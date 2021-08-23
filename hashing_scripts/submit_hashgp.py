#script to submit a job for each test data-set
import os
import sys
import time


files = [
#  "lsh_testdata_2cl"
#  "cellmix_rnaseq_traj",
#  "rnamixcs2_rnaseq_traj",
#  "rnamixss_rnaseq_traj",
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl",
  "cellbench10x_rnaseq_3cl",
  "cellbench10x_rnaseq_5cl",
  "cellbenchDS_rnaseq_3cl",
  "cellbenchCS2_rnaseq_3cl",
  "cellbenchCS2_rnaseq_5cl"
]

gsets = [
  "_full",
  "_hvg5k",
  "_hvg1k",
  "_hvg500"
]

fnames = [x+y for x in files for y in gsets]
idx_test = ['cl1','jsd','l1','l2','s2jsd']

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/calc_hashgp.R '
data_dir  = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/tables/hashes/ '


for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i)
    name_tmp =  x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.calcgp.err '
    out = '-o ' + out_dir + name_tmp + '.calcgp.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    #print(args)
    os.system(cmd)
    time.sleep(1)
