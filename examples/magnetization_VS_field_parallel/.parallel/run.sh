#!/bin/bash
#SBATCH -n 1
#SBATCH -t 1:1:00
#SBATCH --mem-per-cpu=5000
#SBATCH --exclude=milan[8,10,23]
#SBATCH --array=0-39
srun python run.py
