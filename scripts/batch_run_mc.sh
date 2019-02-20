#!/bin/bash -l
#SBATCH --job-name=run_mc
#SBATCH --nodes=1
#SBATCH --constraint=mc
#SBATCH --time=01:30:00
#SBATCH --output=run_mc.out
#SBATCH --error=run_mc.err

./run-bench.sh arbor --model "ring kway" --config "small large"
