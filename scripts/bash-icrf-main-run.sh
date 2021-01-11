for i in {1..6}; do for j in {1..300}; do sbatch --time=5:00:00 bash-icrf-sim.sh $i $j 1 0; done; done
for i in {1..6}; do for j in {1..300}; do sbatch --time=5:00:00 bash-icrf-sim.sh $i $j 3 0; done; done
# for j in {1..100}; do sbatch bash-icrf-sim.sh 4 $j 1 0; done
# for j in {1..100}; do sbatch bash-icrf-sim.sh 4 $j 3 0; done
# i: scenario, j: sim replicates, k: n.monitor, l: 1 for pilot 0 for real.