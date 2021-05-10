#!/bin/bash
#SBATCH --time=30:00:00 --mem=10000
module load r/3.5.2
for i in {1..6}; 
  do Rscript --vanilla scripts/1A.main.R $i $1 $2 0 $3 #$1 = j (sim), $2 = n.monitor, $3 = date
done