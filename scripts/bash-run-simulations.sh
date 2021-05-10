declare -a date=(2021-05-10)

for j in {1..300}; do 
  echo "1A. bash-icrf-main.sh"
  sbatch --time=30:00:00 scripts/bash-icrf-main.sh $j 1 $date
  sbatch --time=30:00:00 scripts/bash-icrf-main.sh $j 3 $date
  
  echo "1B. bash-icrf-sizeSim.sh"
  for i in 100 200 400; do sbatch --mem=10000 --time=1:00:00 scripts/bash-icrf-largeSample.sh $i $j 1 0 $date; done
  # For 1600 a larger memory is needed.
  for i in 800 1600; do sbatch --mem=30000 --time=2:00:00 scripts/bash-icrf-largeSample.sh $i $j 1 0 $date; done
  # i: ntrain, j: sim replicates, k: n.monitor, l: 1 for pilot 0 for real.
  
  echo "1C. bash-icrf-honesty.sh"
  sbatch scripts/bash-icrf-honesty.sh $j $date
  
  echo "5. bash-avalanche.sh"
  sbatch scripts/bash-avalanche.sh $j
  
  echo "6. bash-nlms.sh"
  sbatch scripts/bash-nlms.sh $j
done