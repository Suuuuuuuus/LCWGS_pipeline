#!/bin/bash
#SBATCH -J perform_genotype_calling_on_1x_samples
#SBATCH -A ansari.prj
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/ansari/users/gjx698/low_coverage_sequencing/scripts
#SBATCH -c 4
#SBATCH --array=1-22

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

zgrep -v '#' "../data/gnomAD_vcf/gnomad.genomes.v3.1.2.sites.chr${SLURM_ARRAY_TASK_ID}.vcf.bgz" | \
awk -v "OFS=\t" '{ $8=$8; sub(/^.*AF=/, "", $8); sub(/;.*/, "", $8); print $1,$2,$4,$5,$8 }' | \
awk '$5 !~ /^[[:alpha:]]/' | \
awk '!(length($3)>1 || length($4)>1)' | \
awk '!($5==0.00000)' > "../results/variant_calling/gnomAD_MAFs/gnomAD_MAF_chr${SLURM_ARRAY_TASK_ID}.txt"
