output_date=2024_05_08
parent_dir=/well/band/users/rbx225/GAMCC/data/hla_ref_panel/new/
inputs_dir=${parent_dir}quilt_hla_${output_date}/
test_dir=${parent_dir}quilt_hla_${output_date}/
reference_package_dir=${inputs_dir}quilt_hla_reference_panel_files/
mkdir -p ${test_dir}
mkdir -p ${inputs_dir}
mkdir -p ${reference_package_dir}

cd ${inputs_dir}
oneKG_vcf_name=CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.filtered.shapeit2-duohmm-phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${oneKG_vcf_name}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${oneKG_vcf_name}.tbi

bcftools view --output-file oneKG.temp.vcf.gz --output-type z --min-alleles 2 --max-alleles 2 --types snps ${oneKG_vcf_name} chr6:25000000-34000000
tabix oneKG.temp.vcf.gz
bcftools convert --haplegendsample oneKG oneKG.temp.vcf.gz
reference_haplotype_file=${inputs_dir}oneKG.hap.gz
reference_legend_file=${inputs_dir}oneKG.legend.gz
reference_sample_file=${inputs_dir}oneKG.samples
## convert to uppercase, slight wording change
cd ..
sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' ${reference_sample_file}

## To download IPD-IGMT version 3.39, for example
ipdigmt_link=https://github.com/ANHIG/IMGTHLA/blob/032815608e6312b595b4aaf9904d5b4c189dd6dc/Alignments_Rel_3390.zip?raw=true

cd ${test_dir}
wget ${ipdigmt_link}
ipdigmt_filename_extra=`basename ${ipdigmt_link}`
ipdigmt_filename=`basename ${ipdigmt_link} | sed 's/?raw=true//g'`
mv ${ipdigmt_filename_extra} ${ipdigmt_filename}

cd ${test_dir}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

cd ${inputs_dir}
wget ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20130507_omni_recombination_rates/ACB_omni_recombination_20130507.tar
tar -xvf ACB_omni_recombination_20130507.tar
R -f /well/band/users/rbx225/software/QUILT/scripts/make_b38_recomb_map.R --args . ACB 6

refseq_table_file=/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz

cd ${inputs_dir}
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
ref_fasta=${inputs_dir}GRCh38_full_analysis_set_plus_decoy_hla.fa

cd /well/band/users/rbx225/software/QUILT/
./QUILT_HLA_prepare_reference.R \
--outputdir=${reference_package_dir} \
--nGen=100 \
--hla_types_panel=${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt \
--ipd_igmt_alignments_zip_file=${test_dir}${ipdigmt_filename} \
--ref_fasta=${ref_fasta} \
--refseq_table_file=${refseq_table_file} \
--full_regionStart=25587319 \
--full_regionEnd=33629686 \
--buffer=500000 \
--region_exclude_file=hla_ancillary_files/hlagenes.txt \
--genetic_map_file=${inputs_dir}ACB/ACB-chr6-final.b38.txt.gz \
--reference_haplotype_file=${reference_haplotype_file} \
--reference_legend_file=${reference_legend_file} \
--reference_sample_file=${reference_sample_file} \
--reference_exclude_samples_for_initial_phasing=FALSE \
--hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
--nCores=6
