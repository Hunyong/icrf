declare -a arr=(100 200 400 800)
for i in "${arr[@]}"; do for j in {1..300}; do sbatch --mem=10000 --time=6:00:00 bash-icrf-sizeSim.sh $i $j 1 0; done; done

# For 1600 a larger memory is needed.
for j in {1..300}; do sbatch --mem 30000 --time 45:00:00 bash-icrf-sizeSim.sh 1600 $j 1 0; done
# i: ntrain, j: sim replicates, k: n.monitor, l: 1 for pilot 0 for real.