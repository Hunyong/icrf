#!/bin/bash
#SBATCH --time=3:00:00 --mem=10000
module load r/3.5.2
Rscript --vanilla ./5avalanche.R $1 $2  #1: i (replicate), 2: nfold (=10 by default)