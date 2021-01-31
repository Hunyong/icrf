#!/bin/bash
#SBATCH --time=3:00:00 --mem=20000
module load r/3.5.2
Rscript --vanilla scripts/5avalanche-noise.R $1 $2  #1: i (replicate: 1 ~ 300), 2: nfold (=10 by default)