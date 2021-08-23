#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -cwd
#$ -e ./outputs/lunglymph.err
#$ -o ./outputs/lunglymph.out
#$ -l h_fsize=20G   
module load conda_R

Rscript ./scripts/lunglymph_PoiDist.R 
