#!/bin/bash

#SBATCH -J $JOBNAME
#SBATCH --partition DPB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=6144
#SBATCH --share
#SBATCH --time=1000:00:00
#SBATCH --mail-type=FAIL

source /etc/profile
ulimit -t unlimited

srun $CMD
