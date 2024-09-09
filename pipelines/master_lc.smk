include: "preprocess.smk"
#include: "fastqc.smk"
#include: "reference.smk"
#include: "alignment.smk"

include: "subsample.smk"
#include: "kmer.smk"
include: "dup_rate.smk"
include: "coverage.smk"

include: "imputation_prep.smk"
include: "imputation.smk"

include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

concatenate = config['concatenate']
RECOMB_POP = config["RECOMB_POP"]
panels = config['panels']

# panels are a list of panels, whereas PANEL_NAME is the current panel in use
# decide later whether to incorporate into the pipeline

rule preprocess_all:
    input:
        fwd_pair = expand("data/fastq_cleaned/{id}_1.fastq.gz", id = samples_lc),
        rev_pair = expand("data/fastq_cleaned/{id}_2.fastq.gz", id = samples_lc),
        fwd_unpair = expand("data/fastq_cleaned/{id}_unpaired_1.fastq.gz", id = samples_lc),
        rev_unpair = expand("data/fastq_cleaned/{id}_unpaired_2.fastq.gz", id = samples_lc)

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "data/ref_panel/malariaGen_v3_b38/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

rule reference_all:
    input:
        amb = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.amb" if concatenate else "data/references/GRCh38.fa.amb",
        ann = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.ann" if concatenate else "data/references/GRCh38.fa.ann",
        bwt = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.bwt" if concatenate else "data/references/GRCh38.fa.bwt",
        pac = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.pac" if concatenate else "data/references/GRCh38.fa.pac",
        sa = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.sa" if concatenate else "data/references/GRCh38.fa.sa",

        v3_b37 = expand("data/ref_panel/malariaGen_v3_b37_alone/malariaGen_v3_b37_alone.chr{chr}.vcf.gz", chr = chromosome),
        v1_b38 = expand("data/ref_panel/malariaGen_v1_b38/malariaGen_v1_b38.chr{chr}.vcf.gz", chr = chromosome),
        v3_b38_alone = expand("data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz", chr = chromosome),
        dictionary = "data/references/GRCh38_with_alt.dict",

        mGen_chunk_vcfs = [mGen_chunk_vcf_lst],
        v3_b38 = expand("data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz", chr = chromosome)

rule fastqc_all:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html", id = samples_lc),
        html2 = expand("results/fastqc/{id}_2_fastqc.html", id = samples_lc),
        zip1 = expand("results/fastqc/{id}_1_fastqc.zip", id = samples_lc),
        zip2 = expand("results/fastqc/{id}_2_fastqc.zip", id = samples_lc),
        multiqc_lc = "results/fastqc/multiqc_lc/multiqc_report.html",

        error_rate_tsv = "results/fastqc/per_base_error_rate.tsv",
        bqsr_reports = expand("results/fastqc/BQSR_reports/{id}.BQSR.report", id = samples_lc)

rule alignment_all:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_lc),
        bais = expand("data/bams/{id}.bam.bai", id = samples_lc)

rule subsample_all:
    input:
        ss_bams = expand("data/subsampled_bams/{id}_subsampled.bam", id = samples_lc),
        ss_bais = expand("data/subsampled_bams/{id}_subsampled.bam.bai", id = samples_lc),
        metrics = expand("data/bams/{id}.metrics.txt", id = samples_lc)

rule coverage_all:
    input:
        # bedgraphs = expand("results/coverage/bedgraphs/{id}_bedgraph_nozero.bed", id = samples_lc),
    #    ss_bedgraphs = expand("results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.bed", id = samples_lc),
    #    cumsum_ary = expand("results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt", id = samples_lc),
    #    per_bin_coverage_1x_coordinates = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt", id = samples_lc, chr = chromosome),
    #    per_bin_coverage_1x_bases = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt", id = samples_lc, chr = chromosome),
    #    per_chromosome_coverage = expand("results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt", id = samples_lc),
    #    ss_per_chromosome_coverage = expand("results/coverage/per_chromosome_ss_coverage/{id}_per_chromosome_ss_coverage.txt", id = samples_lc),
    #    uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
    #    ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt"
        # avg_coverage = "results/coverage/per_sample_coverage.txt",
        hc_bedgraphs = expand("results/coverage/hc_bedgraphs/{hc}_bedgraph_nozero.bed", hc = samples_hc),
        hc_coverage = "results/coverage/hc_coverage.txt"

rule dup_rate_all:
    input:
        samtools = "results/dup_rate/duplication_rate_samtools.txt",
        avg_fragment_size = "results/fragment_size/fragment_size.txt",
        proportion_ss_fragment_size = "results/fragment_size/proportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/proportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"

rule kmer_all:
    input:
        jf_read = expand("results/kmer/{id}/read{read}/{id}_read{read}.tsv", id = samples_lc, read = ['1','2']),
        jf_quality = expand("results/kmer/{id}/read{read}/{id}_quality{read}.tsv", id = samples_lc, read = ['1','2']),
        jf_position = expand("results/kmer/{id}/read{read}/{id}_position{read}.tsv", id = samples_lc, read = ['1','2']),
        kmer_accuarcy1 = "results/kmer/kmer_accuracy_read1.txt",
        kmer_accuarcy2 = "results/kmer/kmer_accuracy_read2.txt",
        per_bin_kmer_accuracy1 = "results/kmer/per_bin_kmer_accuracy_read1.txt",
        per_bin_kmer_accuracy2 = "results/kmer/per_bin_kmer_accuracy_read2.txt",
        fragment_length = expand("results/kmer/{id}/fragment_length.tsv", id = samples_lc),
        per_bin_kmer_accuracy = expand("results/kmer/{id}/read{read}/per_bin_kmer_error_rate_read{read}.txt", id = samples_lc, read = ['1', '2']),
        graph_kmer_position = expand("graphs/kmer_position/{id}_kmer_position.png", id = samples_lc)

rule imputation_prep_all:
    input:
        bamlist = "results/imputation/bamlist.txt",
        recomb = expand("results/imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        json = expand("results/imputation/refs/{panel}/regions.json", panel = panels),
        hap = expand("results/imputation/refs/{panel}/{panel}.chr{chr}.hap.gz", chr = chromosome, panel = panels),
        legend = expand("results/imputation/refs/{panel}/{panel}.chr{chr}.legend.gz", chr = chromosome, panel = panels),
        samples = expand("results/imputation/refs/{panel}/{panel}.chr{chr}.samples", chr = chromosome, panel = panels)

all_RData = {}
all_vcf_lst = {}
all_vcf_dict = {}
for p in panels:
    region = "results/imputation/refs/" + p + "/regions.json"
    ref_prefix = "results/imputation/refs/" + p + "/RData/ref_package.chr"
    vcf_prefix = "results/imputation/vcfs/" + p + "/regions/quilt.chr"
    all_RData[p], all_vcf_lst[p], all_vcf_dict[p] = get_vcf_concat_lst(region, ref_prefix, vcf_prefix)

rule imputation_all:
    input:
        RData = [convert_dict_to_lst(all_RData)],
        vcf_regions = [convert_dict_to_lst(all_vcf_lst)],
        vcfs = expand("results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz", panel = panels, chr = chromosome)
