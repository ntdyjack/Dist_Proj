#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_fsize=20G
#$ -pe local 1
#$ -q shared.q@compute-08[5-9],shared.q@compute-09[0-9],shared.q@compute-10[0-9],shared.q@compute-11[0-6]
ml conda_R
Rscript $@
