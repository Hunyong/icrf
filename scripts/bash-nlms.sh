#!/bin/bash
#SBATCH --time=1:00:00 --mem=30000
module load r/3.5.2
Rscript --vanilla scripts/6NLMS.R $1 #1: i (replicate: 1 ~ 300)