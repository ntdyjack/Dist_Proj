#script to submit a job for each test data-set
import os
import sys
import time

files = [
  #"pbmc_citeseq_5cl"
  "pbmc_citeseq_5cl_hvg500",
  "pbmc_citeseq_5cl_hvg1k",
  #"malt_citeseq_2cl",
  "malt_citeseq_2cl_hvg500",
  "malt_citeseq_2cl_hvg1k",
  #"jurkat293t_rnaseq_2cl_"
  "jurkat293t_rnaseq_2cl_hvg500",
  "jurkat293t_rnaseq_2cl_hvg500",
]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_rscript.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/disttests.R '
data_dir  = '/users/ndyjack/Dist_Proj/tables/ '
output_dir = '/users/ndyjack/Dist_Proj/outputs/'

for x in files:
    args = sh_script + r_script + data_dir + output_dir + ' ' + x
    err = '-e ' + output_dir + x + '.disttest.err '
    out = '-o ' + output_dir + x + '.disttest.out '
    cmd = 'qsub -N ' + x + ' ' + err + out + ' ' + args
#    print(cmd)
    os.system(cmd)
    time.sleep(2)
