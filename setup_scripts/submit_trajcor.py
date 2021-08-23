#script to submit a job for each test data-set
import os
import sys
import time


#files = [
#  "cellmix_rnaseq_traj",
#  "rnamixcs2_rnaseq_traj",
#  "rnamixss_rnaseq_traj"
#]

#gsets = [
#  "_full",
#  "_hvg5k",
#  "_hvg1k",
#  "_hvg500"
#]

#fnames = [x+y for x in files for y in gsets]

fnames = os.listdir('/fastscratch/myscratch/ndyjack/simulations/')
fnames = [t for t in fnames if '_t' in t]

idx_test = [0,1,2,3,4,5,6,7,8,9,10,12,13,14,23]

sh_script = '/users/ndyjack/Dist_Proj/scripts/setup_scripts/qsub_gplus.sh '
#r_script = '/users/ndyjack/Dist_Proj/scripts/figure4_scripts/calc_trajcor.R '
r_script = '/users/ndyjack/Dist_Proj/scripts/figure4_scripts/calc_trajcor_splat.R '
#data_dir  = '/users/ndyjack/Dist_Proj/pickles/ '
data_dir  = '/fastscratch/myscratch/ndyjack/pickles/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '

for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i) + ' '
    name_tmp = x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.trajcor.err '
    out = '-o ' + out_dir + name_tmp + '.trajcor.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    #print(args)
    os.system(cmd)
    time.sleep(1)
