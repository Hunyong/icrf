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
