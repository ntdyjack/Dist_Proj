#!/bin/bash
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_fsize=10G
#$ -pe local 1
#$ -q shared.q@compute-08[5-9],shared.q@compute-09[0-9],shared.q@compute-10[0-9],shared.q@compute-11[0-6]
source activate clean
python3 $@
