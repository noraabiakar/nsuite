#!/bin/bash -l
#SBATCH --job-name=build_test
#SBATCH --nodes=1
#SBATCH --constraint=mc
#SBATCH --time=00:15:00
#SBATCH --output=build.out

./install-local.sh arbor -e ./systems/daint-mc.sh

