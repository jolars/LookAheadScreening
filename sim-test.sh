#!/bin/sh

#SBATCH -t 00:15:00

#SBATCH --mail-user=johan.larsson@stat.lu.se
#SBATCH --mail-type=ALL

#SBATCH -J test
#SBATCH -o test_%j.out
#SBATCH -e test_%j.err

#SBATCH -N 1
#SBATCH --tasks-per-node=20

# modules
module purge

# bind results folders in Project to host
export SINGULARITY_BIND="results:/Project/results"

# run the test
singularity run container.sif test.R