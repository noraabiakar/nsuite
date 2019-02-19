#!/bin/bash -l
#SBATCH --job-name=run_test
#SBATCH --nodes=1
#SBATCH --constraint=mc
#SBATCH --time=00:15:00
#SBATCH --output=run.out

./run-bench.sh arbor --model "ring" --config "small"
