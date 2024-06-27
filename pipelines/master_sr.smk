include: "short_read.smk"

#include: "test.smk"
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

chromosome = [i for i in range(1,23)]

read_lengths = ['151-long', '151-optimal', '151-real', '151-short', '300-long', '300-optimal', '300-real', '300-short']
means = ['500', '500', '329', '329', '1000', '1000', '830', '830']
stds = ['13', '13', '13', '13', '200', '200', '200', '200']
error1 = ['0.0024-0.0071', '0', '0.0024-0.0071', '0', '0.0001', '0', '0.0001', '0']
error2 = ['0.0034-0.0105', '0', '0.0034-0.0105', '0', '0.0001', '0', '0.0001', '0']
haplotypes = ['mat', 'pat']
coverage = '0.6'

RECOMB_POP = config["RECOMB_POP"]
PANEL_NAME = config["hc_panel"]

rule short_read_all:
    input:
        fastq1 = expand("data/sr_simulations/{rl}/{rl}.fastq1.gz", rl = read_lengths),
        fastq2 = expand("data/sr_simulations/{rl}/{rl}.fastq2.gz", rl = read_lengths),
        bams = expand("data/sr_bams/{rl}.bam", rl = read_lengths),
        bais = expand("data/sr_bams/{rl}.bam.bai", rl = read_lengths)

region_file = "data/imputation_accessories/5Mb_chunks.json"
region = "results/sr_imputation/refs/" + PANEL_NAME + "/regions.json"
ref_prefix = "results/sr_imputation/refs/" + PANEL_NAME + "/RData/ref_package.chr"
vcf_prefix = "results/sr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr"

sr_oneKG_RData, sr_oneKG_vcf_lst, sr_oneKG_vcf_dict = get_vcf_concat_lst(region, ref_prefix, vcf_prefix)

rule sr_imputation_prep_all:
    input:
        bamlist = "results/sr_imputation/bamlist.txt",
        hap = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        legend = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        samples = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

rule sr_imputation_all:
    input:
        RData = [sr_oneKG_RData],
        vcf_regions = [sr_oneKG_vcf_lst],
        vcfs = expand(f"results/sr_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz", chr = chromosome),
        truth = expand("results/sr_imputation/truth/short_read_truth.chr{chr}.vcf.gz", chr = chromosome)

pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['sr_imputation_dir']
lc_vcf_dir = config['sr_lc_vcf_dir']
hc_vcf_dir = config['sr_hc_vcf_dir']

rule sr_imputation_comparison_all:
    input:
        h_report_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_report_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_impacc_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc = expand("{imp_dir}impacc/all_samples/by_variant/all.h.impacc.tsv", imp_dir = imputation_dir),
        v_impacc = expand("{imp_dir}impacc/all_samples/by_sample/all.v.impacc.tsv", imp_dir = imputation_dir),

        lc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        hc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        afs = expand("{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv", imp_dir = imputation_dir, chr = chromosome),
        
        sumstats = expand("{imp_dir}summary_metrics.tsv", imp_dir = imputation_dir)