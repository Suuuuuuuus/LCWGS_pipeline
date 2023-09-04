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

id=$(cat ../samples.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

dir1="../../GAMCC_1x/data/fastq_cleaned/WTCHG_991075_"
dir2="../../GAMCC_0.1x/data/fastq_cleaned/WTCHG_988637_"

mkdir -p ../data/merged_fastqs
zcat "$dir1$id_1.fastq.gz" > "../data/merged_fastqs/$id_1.fastq"
zcat "$dir2$id_1.fastq.gz" >> "../data/merged_fastqs/$id_1.fastq"
gzip "../data/merged_fastqs/$id_1.fastq"
zcat "$dir1$id_2.fastq.gz" > "../data/merged_fastqs/$id_2.fastq"
zcat "$dir2$id_2.fastq.gz" >> "../data/merged_fastqs/$id_2.fastq"
gzip "../data/merged_fastqs/$id_2.fastq"
