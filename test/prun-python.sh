#!/bin/bash --login

#SBATCH --job-name=porphyrin
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00

export OPENBLAS_NUM_THREADS=1
#export OMP_NUM_THREADS=36
#export LD_LIBRARY_PATH=/home/vposligua/miniconda3/lib/
#export PATH=/home/vposligua/bin/xtb/compiled/usr/local/bin:${PATH}

ulimit -s unlimited
srun python3 polymer.py
