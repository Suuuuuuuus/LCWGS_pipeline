#!/bin/bash
#SBATCH -J perform_genotype_calling_on_1x_samples
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/band/users/rbx225/recyclable_files
#SBATCH -c 4
#SBATCH --array=1-22

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

mkdir -p "AFs/oneKG_MAFs/"

zgrep -v '#' "ref_panels/oneKG/oneKG.chr${SLURM_ARRAY_TASK_ID}.vcf.gz" | \
awk -v "OFS=\t" '{ $8=$8; sub(/^.*AF_AFR=/, "", $8); sub(/;.*/, "", $8); print $1,$2,$4,$5,$8 }' | \
#awk '$5 !~ /^[[:alpha:]]/' | \
#awk '!(length($3)>1 || length($4)>1)' | \
#awk '!($5==0.00000)' > "results/variant_calling/gnomAD_MAFs/gnomAD_MAF_chr${SLURM_ARRAY_TASK_ID}.txt"
awk '$5 !~ /^[[:alpha:]]/' > "AFs/oneKG_MAFs/oneKG_AF_AFR_chr${SLURM_ARRAY_TASK_ID}.txt"
