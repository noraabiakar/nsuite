#!/bin/bash -l
#SBATCH --job-name=run_gpu
#SBATCH --nodes=1
#SBATCH --constraint=gpu
#SBATCH --time=01:30:00
#SBATCH --output=run_gpu.out
#SBATCH --error=run_gpu.err

./run-bench.sh arbor --model "ring kway" --config "small large"
