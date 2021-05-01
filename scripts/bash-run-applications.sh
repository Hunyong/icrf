echo "bash-avalanche.sh";
for j in {1..300}; do sbatch scripts/bash-avalanche.sh $j; done
echo "bash-nlms.sh";
for j in {1..300}; do sbatch scripts/bash-nlms.sh $j; done
