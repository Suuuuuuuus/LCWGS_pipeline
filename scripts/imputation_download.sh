#!/bin/bash
##SBATCH -J download_files_for_imputation
##SBATCH -A ansari.prj
##SBATCH -p short
##SBATCH -o output.out
##SBATCH -D /well/band/users/rbx225/GGVP/
##SBATCH -c 4
##SBATCH --array=1-22

#### Download relevant data for imputation:
# These data should be downloaded manually if not supplied. The cluster the developer uses has no internet access on its nodes.
# One should run this script in the top-level folder

RECOMB_POP=ACB # Note that GWD (Gambian) is not available from the ncbi website

cd results/imputation/
wget "ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/${RECOMB_POP}_omni_recombination_20130507.tar"
tar -xvf "${RECOMB_POP}_omni_recombination_20130507.tar"
rm "${RECOMB_POP}_omni_recombination_20130507.tar"

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz