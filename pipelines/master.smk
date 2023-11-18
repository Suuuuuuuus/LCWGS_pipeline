#include: "preprocess.smk"
#include: "alignment.smk"
#include: "kmer.smk"
#include: "dup_rate.smk"
#include: "rmdup.smk"
#include: "fastqc.smk"
#include: "subsample.smk"
#include: "reference.smk"
include: "coverage.smk"
include: "imputation.smk"
#include: "imputation_prep.smk"
#include: "variant_calling.smk"

configfile: "pipelines/config.json"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
sys.path.append("scripts")
import lcwgSus

# samples = pd.read_table(config['samples'], header = None, names = ['Code'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)
panels = config['panels']

chromosome = [i for i in range(1,23)]

# The followings are global parameters:
clean_fastq = config['clean_fastq']
reheader = config['reheader']
concatenate = config['concatenate']

# The followings are global parameters from `activate`:
QUILT_HOME = config["QUILT_HOME"]
ANALYSIS_DIR = config["ANALYSIS_DIR"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
NGEN=config["NGEN"]
RECOMB_POP=config["RECOMB_POP"]
PANEL_NAME=config["PANEL_NAME"]

rule variant_calling_all:
    input:
        imputation_vcf = "results/imputation/tmp/res.txt"

rule preprocess_all:
    input:
        fwd_pair = expand("data/fastq_cleaned/{id}_1.fastq.gz", id = ids_1x_all),
        rev_pair = expand("data/fastq_cleaned/{id}_2.fastq.gz", id = ids_1x_all),
        fwd_unpair = expand("data/fastq_cleaned/{id}_unpaired_1.fastq.gz", id = ids_1x_all),
        rev_unpair = expand("data/fastq_cleaned/{id}_unpaired_2.fastq.gz", id = ids_1x_all)

rule alignment_all:
    input:
        bams = expand("data/bams/{id}.bam", id = ids_1x_all),
        bais = expand("data/bams/{id}.bam.bai", id = ids_1x_all)

rule reference_all:
    input:
        amb = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.amb" if concatenate else "data/references/GRCh38.fa.amb",
        ann = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.ann" if concatenate else "data/references/GRCh38.fa.ann",
        bwt = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.bwt" if concatenate else "data/references/GRCh38.fa.bwt",
        pac = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.pac" if concatenate else "data/references/GRCh38.fa.pac",
        sa = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.sa" if concatenate else "data/references/GRCh38.fa.sa"

rule subsample_all:
    input:
        ss_bams = expand("data/subsampled_bams/{id}_subsampled.bam", id = ids_1x_all),
        ss_bais = expand("data/subsampled_bams/{id}_subsampled.bam.bai", id = ids_1x_all)

rule coverage_all:
    input:
#        bedgraphs = expand("results/coverage/bedgraphs/{id}_bedgraph.bed", id = ids_1x_all),
#        ss_bedgraphs = expand("results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.bed", id = ids_1x_all),
#        cumsum_ary = expand("results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt", id = ids_1x_all),
#        per_bin_coverage_1x_coordinates = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt", id = ids_1x_all, chr = chromosome),
#        per_bin_coverage_1x_bases = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt", id = ids_1x_all, chr = chromosome),
#        per_chromosome_coverage = expand("results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt", id = ids_1x_all),
#        ss_per_chromosome_coverage = expand("results/coverage/per_chromosome_ss_coverage/{id}_per_chromosome_ss_coverage.txt", id = ids_1x_all),
#        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
        ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt"
#        avg_coverage = "results/coverage/per_sample_coverage.txt"

rule dup_rate_all:
    input:
        samtools = "results/dup_rate/duplication_rate_samtools.txt",
        avg_fragment_size = "results/fragment_size/fragment_size.txt",
        proportion_ss_fragment_size = "results/fragment_size/proportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/proportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"

rule rmdup_all:
    input:
        dedup_bams = expand("data/dedup_bams/{id}.bam", id = ids_1x_all),
        dedup_bais = expand("data/dedup_bams/{id}.bam.bai", id = ids_1x_all)

rule fastqc_all:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html", id = ids_1x_all),
        html2 = expand("results/fastqc/{id}_2_fastqc.html", id = ids_1x_all),
        zip1 = expand("results/fastqc/{id}_1_fastqc.zip", id = ids_1x_all),
        zip2 = expand("results/fastqc/{id}_2_fastqc.zip", id = ids_1x_all),
        fastqc = "results/fastqc/duplication_rate_fastqc.txt",
        multiqc = "results/fastqc/multiqc_report.html",
        multiqcdir = "results/fastqc/multiqc_data"

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

# Dumps
REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if exists(file):
    print("Replacing regions to impute with derived file")
    with open(file) as json_file:
        REGIONS = json.load(json_file)

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
        file="results/imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)

rule imputation_prep_all:
    input:
        bamlist = "results/imputation/bamlist.txt",
        recomb = expand("results/imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        json = "results/imputation/regions.json",
        hap = expand(f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        legend = expand(f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        samples = expand(f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

chip_to_extract = pd.read_table("chip.tsv", header = None).loc[:, 0].to_list()
seq_to_extract = sample_linker[sample_linker['Chip_Name'].isin(chip_to_extract)].Seq_Name.to_list()

rule imputation_all:
    input:
        RData = [regions_to_prep],
        vcf_regions = [vcfs_to_impute],
        vcfs = [final_vcfs],
        r2 = expand("results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv", id = seq_to_extract, panel = panels)
        #graph = "graphs/imputation_vs_chip.png"

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
        proportion_ss_fragment_size = "results/fragment_size/proportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/proportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"
    output:
        result = "results/lcwgs_results.csv"
    shell: """
        python {input.script}
    """

