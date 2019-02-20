#!/bin/bash -l
#SBATCH --job-name=run_tave
#SBATCH --nodes=1
#SBATCH --constraint=flat,quad
#SBATCH --time=01:30:00
#SBATCH --output=run_tave.out
#SBATCH --error=run_tave.err

./run-bench.sh arbor --model "ring kway" --config "small"
