#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -cwd
#$ -l h_fsize=20G
#$ -pe local 3
module load conda_R
#Rscript /users/ndyjack/Dist_Proj/scripts/splatterTests_v2.R $1
#Rscript /users/ndyjack/Dist_Proj/scripts/splatterTests_dr.R
Rscript /users/ndyjack/Dist_Proj/scripts/splatterTests_nG.R $1
## -e /users/ndyjack/Dist_Proj/outputs/splatterTests.nG.err
## -o /users/ndyjack/Dist_Proj/outputs/splatterTests.nG.out
