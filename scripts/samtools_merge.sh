#!/bin/bash
#SBATCH -J samtools_merge
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/band/users/rbx225/GAMCC_1x/scripts
#SBATCH -c 8
#SBATCH --array=1-2

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

full=$(cat ../samples.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
id="${full##*_}"

dir1="../data/bams/WTCHG_991075_"
dir2="../../GAMCC_0.1x/data/bams/WTCHG_988637_"

mkdir -p ../data/merged_bams
samtools merge -o "../data/merged_bams/$id.bam" "$dir1$id.bam" "$dir2$id.bam"
samtools index "../data/merged_bams/$id.bam"
