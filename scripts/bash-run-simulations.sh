echo "1. bash-icrf-main.sh"
for j in {1..300}; do sbatch --time=30:00:00 scripts/bash-icrf-main.sh $j 0; done

echo "2. bash-icrf-honesty.sh"
for j in {1..300}; do sbatch scripts/bash-icrf-honesty.sh $j; done

echo "3. bash-icrf-sizeSim.sh"
declare -a arr=(100 200 400)
declare -a arr2=(800 1600)
for i in "${arr[@]}"; do for j in {1..300}; do sbatch --mem=10000 --time=6:00:00 scripts/bash-icrf-sizeSim.sh $i $j 1 0; done; done
# For 1600 a larger memory is needed.
for i in "${arr2[@]}"; do for j in {1..300}; do sbatch --mem=50000 --time=45:00:00 scripts/bash-icrf-sizeSim.sh $i $j 1 0; done; done
# i: ntrain, j: sim replicates, k: n.monitor, l: 1 for pilot 0 for real.


echo "4. bash-icrf-time.sh"
declare -a arr=(100 200 400)
declare -a arr2=(800 1600)
# n ntree nfold sim
for ntree in 100 300; do
  echo "ntree = $ntree";  
  for nfold in 1 3 5 10; do
    echo "ntree = $ntree, nfold = $nfold";  
    for n in "${arr[@]}"; do 
      echo "ntree = $ntree, nfold = $nfold, n = $n";  
      for j in {1..30}; do sbatch --mem=10000 --time=6:00:00 scripts/bash-icrf-time.sh $n $ntree $nfold $j; done; done
    # For 1600 a larger memory is needed.
    for n in "${arr2[@]}"; do 
      echo "ntree = $ntree, nfold = $nfold, n = $n";  
      for j in {1..30}; do sbatch --mem=50000 --time=45:00:00 scripts/bash-icrf-time.sh $n $ntree $nfold $j; done; done
  done
done  
