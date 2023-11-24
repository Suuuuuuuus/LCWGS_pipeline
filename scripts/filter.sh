#!/bin/bash
#SBATCH -J samtools_merge
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/band/users/rbx225/GAMCC/scripts/
#SBATCH -c 8
##SBATCH --array=1-2

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

python filter_chip_afs.py
