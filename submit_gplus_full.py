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
r_script = '/users/ndyjack/Dist_Proj/scripts/gplus_fulldist.R '
data_dir  = '/users/ndyjack/Dist_Proj/tables/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/rdata/ '
out_file = ' gplus_results_July2020.txt'

#idx_test = [ [x for x in range(0,15)], [x for x in range(16,22)], [x for x in range(28,30)]]
#idx_test = [item for sublist in idx_test for item in sublist]
idx_test = [16,17]

for x in files:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i)
    name_tmp = x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp + out_file
    err = '-e ' + out_dir + name_tmp  + '.gplus.full.err '
    out = '-o ' + out_dir + name_tmp + '.gplus.full.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    os.system(cmd)
    time.sleep(2)
