#!/bin/bash
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_fsize=10G
#$ -pe local 1
ml conda_R
Rscript $@
