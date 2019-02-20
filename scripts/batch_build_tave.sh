#!/bin/bash -l
#SBATCH --job-name=build_tave
#SBATCH --nodes=1
#SBATCH --constraint=flat,quad
#SBATCH --time=00:30:00
#SBATCH --output=build_tave.out
#SBATCH --error=build_tave.err

./install-local.sh arbor -e ./systems/tave.sh
