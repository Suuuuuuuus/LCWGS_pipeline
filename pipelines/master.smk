#include: "alignment.smk"
include: "jellyfish.smk"
include: "dup_rate.smk"
include: "fastqc.smk"
include: "subsample.smk"
include: "index_reference.smk"
include: "coverage.smk"
#include: "variant_calling.smk"

configfile: "pipelines/config.json"

#ids_1x = [samples_1x['id'] for samples_1x in config['samples_1x']]
#ids_20x = [samples_20x['id'] for samples_20x in config['samples_20x']]
#ids_NA12878_1x = [NA12878_1x['id'] for NA12878_1x in config['NA12878_1x']]
#ids_NA12878_20x = [NA12878_20x['id'] for NA12878_20x in config['NA12878_20x']]
chromosome = [i for i in range(1,23)]

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

#ids_1x_all = ids_1x + ids_NA12878_1x
#ids_20x_all = ids_20x + ids_NA12878_20x
#ids_all = ids_1x + ids_20x + ids_NA12878_1x + ids_NA12878_20x

rule all:
	input:
		ss_fastq1 = expand("data/subsampled_fastq/{subsample}_subsampled_1.fastq", subsample = ids_1x_all),
		ss_fastq2 = expand("data/subsampled_fastq/{subsample}_subsampled_2.fastq", subsample = ids_1x_all),
        	ss_bams = expand("data/subsampled_bams/{subsample}_subsampled.bam", subsample = ids_1x_all),
        	ss_bais = expand("data/subsampled_bams/{subsample}_subsampled.bam.bai", subsample = ids_1x_all),
        	ss_bedgraphs = expand("results/coverage/subsampled_bedgraphs/{subsample}_subsampled_bedgraph.txt", subsample = ids_1x_all),

        	amb = "data/reference/GRCh38.fa.amb",
        	ann = "data/reference/GRCh38.fa.ann",
        	bwt = "data/reference/GRCh38.fa.bwt",
       	 	pac = "data/reference/GRCh38.fa.pac",
        	sa = "data/reference/GRCh38.fa.sa",

		html1 = expand("results/fastqc/{id}_1_fastqc.html", id = ids_1x_all),
		html2 = expand("results/fastqc/{id}_2_fastqc.html", id = ids_1x_all),
		zip1 = expand("results/fastqc/{id}_1_fastqc.zip", id = ids_1x_all),
		zip2 = expand("results/fastqc/{id}_2_fastqc.zip", id = ids_1x_all),

#		bams = expand("data/bams/{id}.bam", id = ids_1x_all),
#        	bais = expand("data/bams/{id}.bam.bai", id = ids_1x_all),

        	fastqc = "results/fastqc/duplication_rate_fastqc.txt",
        	samtools = "results/dup_rate/duplication_rate_samtools.txt",
        	uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
                avg_coverage = "results/coverage/per_sample_coverage.txt",

        	bedgraphs = expand("results/coverage/bedgraphs/{id}_bedgraph.txt", id = ids_1x_all),

        	per_chromosome_coverage = expand("results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt", id = ids_1x_all),
                per_bin_coverage_1x_coordinates = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt", id = ids_1x_all, chr = chromosome),
        	per_bin_coverage_1x_bases = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt", id = ids_1x_all, chr = chromosome),

		jf_read = expand("results/kmer/{id}/read{read}/{id}_read{read}.tsv.gz", id = ids_1x_all, read = ['1','2']),
        	jf_quality = expand("results/kmer/{id}/read{read}/{id}_quality{read}.tsv.gz", id = ids_1x_all, read = ['1','2']),
        	jf_position = expand("results/kmer/{id}/read{read}/{id}_position{read}.tsv.gz", id = ids_1x_all, read = ['1','2']),

 	        kmer_accuarcy1 = "results/kmer/kmer_accuracy_read1.txt",
 	        kmer_accuarcy2 = "results/kmer/kmer_accuracy_read2.txt",

		lcwgs_wrap_up = "results/lcwgs_results.csv",

#        	graph_subsample_coverage = "graphs/fig8_prop_genome_at_least_coverage.png",
        	graph_samtools_dup_rate = "graphs/samtools_duplication_rate.png",
        	graph_uncoverage_rate = "graphs/uncoverage_rate.png",
        	graph_chromosome_coverage = expand("graphs/fig6_per_bin_coverage_chr{chr}.png", chr = chromosome),
		graph_kmer_position = expand("graphs/kmer_position/{id}_kmer_position.png", id = ids_1x_all)
#		graph = "graphs/NA12878_imputation_accuracy.png",
#		graph_lcwgs = "graphs/lcwgs_imputation_accuracy.png",

rule aggregate_results:
    input:
        script = "scripts/lcwgs_result_wrap_up.py",
        fastqc_dup_rate = "results/fastqc/duplication_rate_fastqc.txt",
        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
        samtools_dup_rate = "results/dup_rate/duplication_rate_samtools.txt",
        kmer_accuracy1 = "results/kmer/kmer_accuracy_read1.txt",
        kmer_accuracy2 = "results/kmer/kmer_accuracy_read2.txt",
        coverage = "results/coverage/per_sample_coverage.txt"
    output:
        result = "results/lcwgs_results.csv"
    shell: """
        python {input.script}
    """

