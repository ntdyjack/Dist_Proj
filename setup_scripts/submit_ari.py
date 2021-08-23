#script to submit a job for each test data-set
import os
import sys
import time


#files = [
#  "pbmc_citeseq_5cl",
#  "malt_citeseq_2cl",
#  "jurkat293t_rnaseq_2cl",
#  "cellbench10x_rnaseq_3cl",
#  "cellbench10x_rnaseq_5cl",
#  "cellbenchDS_rnaseq_3cl",
#  "cellbenchCS2_rnaseq_3cl",
#  "cellbenchCS2_rnaseq_5cl"
#]

#gsets = [
#  "_full",
#  "_hvg5k",
#  "_hvg1k",
#  "_hvg500"
#]

#fnames = [x+y for x in files for y in gsets]

fnames = os.listdir('/fastscratch/myscratch/ndyjack/simulations/')
fnames = [t for t in fnames if not '_t' in t]

#k = '30'
#idx_test = [ [x for x in range(0,22)], [x for x in range(28,30)],["std_pca", "trn_pca", "glm_pca"]]
#idx_test = [item for sublist in idx_test for item in sublist]
#idx_test = [5,6,12,13]
idx_test = [0,1,2,3,4,5,6,7,8,9,10,12,13,14,23]
#idx_test = [5,13,14,23]

sh_script = '/users/ndyjack/Dist_Proj/scripts/qsub_gplus.sh '
r_script = '/users/ndyjack/Dist_Proj/scripts/calc_ari.R '
#data_dir  = '/users/ndyjack/Dist_Proj/pickles/ '
data_dir  = '/fastscratch/myscratch/ndyjack/pickles/ '
out_dir = '/users/ndyjack/Dist_Proj/outputs/'
res_dir = '/users/ndyjack/Dist_Proj/tables/test_datasets/ '

for x in fnames:
  for i in idx_test:
    arg_tmp =  x + ' ' + str(i) + ' '
    name_tmp = x + '.' + str(i) 
    args = sh_script + r_script + data_dir + res_dir + arg_tmp
    err = '-e ' + out_dir + name_tmp  + '.ari.err '
    out = '-o ' + out_dir + name_tmp + '.ari.out '
    cmd = 'qsub -N ' + name_tmp + ' ' + err + out + ' ' + args
    #print(args)
    os.system(cmd)
    time.sleep(1)
