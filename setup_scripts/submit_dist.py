#script to submit a job for each test data-set
import os
import sys
import time


#files = [
  #"splat_sim_1cl",
  #"only293t_rnaseq_1cl",
  #"cellmix_rnaseq_traj",
  #"rnamixcs2_rnaseq_traj",
  #"rnamixss_rnaseq_traj",
  #"pbmc_citeseq_5cl",
  #"malt_citeseq_2cl",
  #"jurkat293t_rnaseq_2cl",
  #"cellbench10x_rnaseq_3cl",
  #"cellbench10x_rnaseq_5cl",
  #"cellbenchDS_rnaseq_3cl",
  #"cellbenchCS2_rnaseq_3cl",
  #"cellbenchCS2_rnaseq_5cl"
#]

#files = ['bg293t_sim_1cl' + "_" + str(x) for x in range(1,31)]
#files = ['splat_bgsim_1cl' + "_" + str(x) for x in range(1,31)]
#files = ['splatsim_1cl']

#gsets = [
#  "_full",
#  "_hvg5k",
#  "_hvg1k",
#  "_hvg500"
#]

#fnames = [x+y for x in files for y in gsets]
fnames = os.listdir('/fastscratch/myscratch/ndyjack/simulations/')


sh_script = '/users/ndyjack/Dist_Proj/scripts/setup_scripts/qsub_mini.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/setup_scripts/calc_dist.py '
#data_dir  = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '
#data_dir  = '/users/ndyjack/Dist_Proj/tables/null_datasets/ '
data_dir = '/fastscratch/myscratch/ndyjack/simulations/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
#res_dir = '/users/ndyjack/Dist_Proj/pickles/ '
res_dir = '/fastscratch/myscratch/ndyjack/pickles/ '

idx_test = [0,1,2,3,4,5,6,7,8,9,10,12,13,14,23]

for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i)
    name_tmp =  x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.calcdist.err '
    out = '-o ' + out_dir + name_tmp + '.calcdist.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    os.system(cmd)
    time.sleep(1)
