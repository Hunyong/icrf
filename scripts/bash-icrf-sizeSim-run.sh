declare -a arr=(100 200 400)
declare -a arr2=(800 1600)
for i in "${arr[@]}"; do for j in {1..300}; do sbatch --mem=10000 --time=6:00:00 scripts/bash-icrf-sizeSim.sh $i $j 1 0; done; done

# For 1600 a larger memory is needed.
for i in "${arr2[@]}"; do for j in {1..300}; do sbatch --mem=50000 --time=45:00:00 scripts/bash-icrf-sizeSim.sh $i $j 1 0; done; done
# i: ntrain, j: sim replicates, k: n.monitor, l: 1 for pilot 0 for real.