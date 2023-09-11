RECOMB_POP="ACB"
cd results/imputation/
wget "ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/${RECOMB_POP}_omni_recombination_20130507.tar"
tar -xvf "${RECOMB_POP}_omni_recombination_20130507.tar"
rm "${RECOMB_POP}_omni_recombination_20130507.tar"

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
