#!/bin/bash -l
#SBATCH --job-name=mc_test
#SBATCH --nodes=1
#SBATCH --constraint=mc
#SBATCH --time=00:15:00

./install-local.sh arbor -e ./systems/daint-mc.sh

