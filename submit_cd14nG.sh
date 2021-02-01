#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -cwd
#$ -l h_fsize=20G
#$ -pe local 5
module load conda_R

#Rscript ./scripts/cd14_TestDist.R
Rscript /users/ndyjack/Dist_Proj/scripts/cd14_nG.R 
