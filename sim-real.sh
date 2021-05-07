#!/bin/sh

#SBATCH -t 8:00:00

#SBATCH --mail-user=johan.larsson@stat.lu.se
#SBATCH --mail-type=ALL

#SBATCH -J realdata
#SBATCH -o realdata_%j.out
#SBATCH -e realdata_%j.err

#SBATCH -N 1
#SBATCH --tasks-per-node=20

# modules
module purge

singularity run --bind results:/Project/results container.sif \
  experiments/realdata.R
