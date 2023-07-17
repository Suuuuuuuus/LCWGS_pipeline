include: "alignment.smk"
include: "kmer.smk"
include: "dup_rate.smk"
include: "fastqc.smk"
include: "subsample.smk"
include: "index_reference.smk"
include: "coverage.smk"
#include: "imputation.smk"
include: "imputation_prep.smk"

configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
#chromosome = [i for i in range(1,23)]
chromosome = [11]

# The followings are global parameters from `activate`:
QUILT_WRAP_HOME = "/well/band/users/rbx225/GGVP/QUILT-wrap/"
QUILT_HOME = "/well/band/users/rbx225/software/QUILT/"
ANALYSIS_DIR = "/well/band/users/rbx225/GGVP/results/imputation/"
RECOMB_POP="ACB"
WINDOWSIZE=5000000
BUFFER=1000000
NGEN=100
RECOMB_POP="ACB"

rule alignment_all:
    input:
        bams = expand("data/bams/{id}.bam", id = ids_1x_all),
        bais = expand("data/bams/{id}.bam.bai", id = ids_1x_all)

rule index_reference_all:
    input:
        amb = "data/reference/GRCh38.fa.amb",
        ann = "data/reference/GRCh38.fa.ann",
        bwt = "data/reference/GRCh38.fa.bwt",
        pac = "data/reference/GRCh38.fa.pac",
        sa = "data/reference/GRCh38.fa.sa"

rule subsample_all:
    input:
        ss_fastq1 = expand("data/subsampled_fastq/{subsample}_subsampled_1.fastq", subsample = ids_1x_all),
        ss_fastq2 = expand("data/subsampled_fastq/{subsample}_subsampled_2.fastq", subsample = ids_1x_all),
        ss_bams = expand("data/subsampled_bams/{subsample}_subsampled.bam", subsample = ids_1x_all),
        ss_bais = expand("data/subsampled_bams/{subsample}_subsampled.bam.bai", subsample = ids_1x_all)

rule coverage_all:
    input:
        bedgraphs = expand("results/coverage/bedgraphs/{id}_bedgraph.txt", id = ids_1x_all),
        ss_bedgraphs = expand("results/coverage/subsampled_bedgraphs/{subsample}_subsampled_bedgraph.txt", subsample = ids_1x_all),
        graph_subsample_coverage = "graphs/fig8_prop_genome_at_least_coverage.png", # This needs investigation
        per_bin_coverage_1x_coordinates = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt", id = ids_1x_all, chr = chromosome),
        per_bin_coverage_1x_bases = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt", id = ids_1x_all, chr = chromosome),
        per_chromosome_coverage = expand("results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt", id = ids_1x_all),
        ss_per_chromosome_coverage = expand("results/coverage/per_chromosome_ss_coverage/{subsample}_per_chromosome_ss_coverage.txt", subsample = ids_1x_all),
        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
        ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt",
        avg_coverage = "results/coverage/per_sample_coverage.txt",
        graph_uncoverage_rate = "graphs/uncoverage_rate.png",
        graph_chromosome_coverage = expand("graphs/fig6_per_bin_coverage_chr{chr}.png", chr = chromosome)

rule dup_rate_all:
    input:
        samtools = "results/dup_rate/duplication_rate_samtools.txt",
        graph_samtools_dup_rate = "graphs/samtools_duplication_rate.png",
        avg_fragment_size = "results/fragment_size/fragment_size.txt",
        proportion_ss_fragment_size = "results/fragment_size/porportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/porportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"

rule fastqc_all:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html", id = ids_1x_all),
	html2 = expand("results/fastqc/{id}_2_fastqc.html", id = ids_1x_all),
	zip1 = expand("results/fastqc/{id}_1_fastqc.zip", id = ids_1x_all),
	zip2 = expand("results/fastqc/{id}_2_fastqc.zip", id = ids_1x_all),
        fastqc = "results/fastqc/duplication_rate_fastqc.txt"

rule kmer_all:
    input:
        jf_read = expand("results/kmer/{id}/read{read}/{id}_read{read}.tsv", id = ids_1x_all, read = ['1','2']),
        jf_quality = expand("results/kmer/{id}/read{read}/{id}_quality{read}.tsv", id = ids_1x_all, read = ['1','2']),
        jf_position = expand("results/kmer/{id}/read{read}/{id}_position{read}.tsv", id = ids_1x_all, read = ['1','2']),
 	kmer_accuarcy1 = "results/kmer/kmer_accuracy_read1.txt",
 	kmer_accuarcy2 = "results/kmer/kmer_accuracy_read2.txt",
        per_bin_kmer_accuracy1 = "results/kmer/per_bin_kmer_accuracy_read1.txt",
	per_bin_kmer_accuracy2 = "results/kmer/per_bin_kmer_accuracy_read2.txt",
	fragment_length = expand("results/kmer/{id}/fragment_length.tsv", id = ids_1x_all),
        per_bin_kmer_accuracy = expand("results/kmer/{id}/read{read}/per_bin_kmer_error_rate_read{read}.txt", id = ids_1x_all, read = ['1', '2']),
        graph_kmer_position = expand("graphs/kmer_position/{id}_kmer_position.png", id = ids_1x_all)

rule imputation_prep_all:
    input:
        bamlist = "results/imputation/bamlist.txt",
        recomb = expand("results/imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        json = "results/imputation/regions.json",
        hap = expand("results/imputation/refs/ggvp.chr{chr}.hap.gz", chr = chromosome),
        legend = expand("results/imputation/refs/ggvp.chr{chr}.legend.gz", chr = chromosome),
        samples = expand("results/imputation/refs/ggvp.chr{chr}.samples", chr = chromosome)

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if exists(file):
    print("Replacing regions to impute with derived file")
    with open(file) as json_file:
        REGIONS = json.load(json_file) ## python is dumb

regions_to_prep=[]
vcfs_to_impute=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/imputation/refs/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file="results/imputation/vcfs/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)

rule imputation_all:
    input:
        RData = [regions_to_prep]

rule all:
	input:
		lcwgs_wrap_up = "results/lcwgs_results.csv"

rule aggregate_results:
    input:
        script = "scripts/lcwgs_result_wrap_up.py",
        fastqc_dup_rate = "results/fastqc/duplication_rate_fastqc.txt",
        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
        samtools_dup_rate = "results/dup_rate/duplication_rate_samtools.txt",
        kmer_accuracy1 = "results/kmer/kmer_accuracy_read1.txt",
        kmer_accuracy2 = "results/kmer/kmer_accuracy_read2.txt",
        coverage = "results/coverage/per_sample_coverage.txt",
        fragment_size = "results/fragment_size/fragment_size.txt",
        proportion_ss_fragment_size = "results/fragment_size/porportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/porportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"
    output:
        result = "results/lcwgs_results.csv"
    shell: """
        python {input.script}
    """

