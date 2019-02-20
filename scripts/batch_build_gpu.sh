#!/bin/bash -l
#SBATCH --job-name=build_gpu
#SBATCH --nodes=1
#SBATCH --constraint=gpu
#SBATCH --time=00:30:00
#SBATCH --output=build_gpu.out
#SBATCH --error=build_gpu.err

./install-local.sh arbor -e ./systems/daint-gpu.sh
