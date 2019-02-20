#!/bin/bash -l
#SBATCH --job-name=build_mc
#SBATCH --nodes=1
#SBATCH --constraint=mc
#SBATCH --time=00:30:00
#SBATCH --output=build_mc.out
#SBATCH --error=build_mc.err

./install-local.sh arbor -e ./systems/daint-mc.sh
