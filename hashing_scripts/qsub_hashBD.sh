#!/bin/bash
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l h_fsize=5G
#$ -pe local 1
#$ -q shared.q@compute-04[3-9],shared.q@compute-05[7-9],shared.q@compute-06[0-9],shared.q@compute-07[2-6]
tmp=$@
rs=${tmp/'/users/ndyjack/Dist_Proj/scripts/hash_10xBD.py'}
ps=${tmp/'/users/ndyjack/Dist_Proj/scripts/hash_10xBD.R'}
ml R/3.6.1 
Rscript $rs
module unload conda/3-4.6.14
module unload conda_R/4.0
python $ps
