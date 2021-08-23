#script to submit a job for each test data-set
import os
import sys
import time


fnames = [
  "full",
  "hvg5k",
  "hvg1k",
  "hvg500"
]

idx_test = [0,1,2,3,4,5,6,7,8,9,10,12,13,14,23]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
#r_script = '/users/ndyjack/Dist_Proj/scripts/calc_gsr.R '
r_script = '/users/ndyjack/Dist_Proj/scripts/calc_gapstats.R '
dt_dir  = '/users/ndyjack/Dist_Proj/pickles/ '
bg_dir = '/fastscratch/myscratch/ndyjack/pickles/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
#res_dir = '/users/ndyjack/Dist_Proj/tables/gsr/ '
res_dir = '/users/ndyjack/Dist_Proj/tables/gapstat/ ' 


for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i) + ' '
    name_tmp = x + '.' + str(i) 
    args = sh_script + r_script + dt_dir + bg_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.gsr.err '
    out = '-o ' + out_dir + name_tmp + '.gsr.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    #print(args)
    os.system(cmd)
    time.sleep(1)
