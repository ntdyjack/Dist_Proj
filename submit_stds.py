#script to submit a job for each test data-set
import os
import sys
import time

files = [
  "pbmc_citeseq_5cl",
  "malt_citeseq_2cl",
  "jurkat293t_rnaseq_2cl"
]

tests = [
  "std_pca",
  "trn_pca",
  "glm_pca"
]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/calc_stds.R '
data_dir  = '/users/ndyjack/Dist_Proj/tables/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/rdata/ '

for x in files:
  for y in tests:
    arg_tmp =  x + ' ' + y
    name_tmp = x + '.' + y
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.gplus.err '
    out = '-o ' + out_dir + name_tmp + '.gplus.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    os.system(cmd)
    time.sleep(2)
