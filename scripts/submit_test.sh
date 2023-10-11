#!/bin/bash
#SBATCH -J samtools_merge
#SBATCH -p short
#SBATCH -o %j.out
#SBATCH -D /well/band/users/rbx225/GAMCC/scripts/
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=1000
##SBATCH --array=1-2

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

python test.py
