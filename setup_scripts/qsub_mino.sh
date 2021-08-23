#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_fsize=20G
#$ -pe local 1
#$ -q shared.q@compute-04[3-9],shared.q@compute-05[7-9],shared.q@compute-06[0-9],shared.q@compute-07[2-6]
python $@
