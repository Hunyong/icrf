#!/bin/bash
#SBATCH --time=3:00:00 --mem=10000
module load r/3.5.2
Rscript --vanilla scripts/1B.largeSample.R $1 $2 $3 $4 $5