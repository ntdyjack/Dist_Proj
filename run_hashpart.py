import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import sys
import numpy as np
import minocore
import h5py


dat = h5py.File('/users/ndyjack/Dist_Proj/tables/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5', 'r')
expr = dat['mm10/data']
