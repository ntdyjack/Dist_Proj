#script to submit a job for each test data-set
import os
import sys
import time

files = [
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl"
]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/gplus_mds.R '
data_dir  = '/users/ndyjack/Dist_Proj/tables/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/rdata/ '

for x in files:
  for i in range(0,17):
    arg_tmp =  x + ' ' + str(i)
    name_tmp = x + '.mds.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.gplus.mds.err '
    out = '-o ' + out_dir + name_tmp + '.gplus.mds.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    os.system(cmd)
    time.sleep(2)
