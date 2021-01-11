for i in {1..6}; do for j in {1..300}; do sbatch bash-icrf-honesty.sh $i $j 1 0; done; done
for i in {1..6}; do for j in {1..300}; do sbatch bash-icrf-honesty.sh $i $j 3 0; done; done