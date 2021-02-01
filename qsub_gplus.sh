#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_fsize=20G
#$ -pe local 1
ml conda_R
Rscript $@
