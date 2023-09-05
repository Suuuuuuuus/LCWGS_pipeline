#!/bin/bash
#SBATCH -J samtools_merge
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/band/users/rbx225/GAMCC/scripts
#SBATCH -c 8
#SBATCH --array=1-241

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

full=$(cat ../samples.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
id="${full##*_}"

mkdir -p ../data/fastq_cleaned
echo 1 > "../data/fastq_cleaned/${id}_1.fastq.gz"
echo 1 > "../data/fastq_cleaned/${id}_unpaired_1.fastq.gz"
echo 1 > "../data/fastq_cleaned/${id}_unpaired_2.fastq.gz"
echo 1 > "../data/fastq_cleaned/${id}_2.fastq.gz"
