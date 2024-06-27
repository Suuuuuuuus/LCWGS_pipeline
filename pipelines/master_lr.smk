include: "long_read.smk"

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

read_lengths = config["lr_read_lengths"]
rls = config["lr_rls"]
haplotypes = config["lr_haplotypes"]
method = config["lr_method"]
coverage = config["lr_sim_coverage"]

RECOMB_POP = config["RECOMB_POP"]
PANEL_NAME = config["hc_panel"]

rule long_read_all:
    input:
        fastqs = expand("data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz", rl = read_lengths, hap = haplotypes),
        bams = expand("data/lr_bams/{rl}.bam", rl = read_lengths),
        bais = expand("data/lr_bams/{rl}.bam.bai", rl = read_lengths),
        fastq = expand("data/lr_simulations/{rl}/{rl}.fastq.gz", rl = read_lengths)

region_file = "data/imputation_accessories/5Mb_chunks.json"
region = "results/lr_imputation/refs/" + PANEL_NAME + "/regions.json"
ref_prefix = "results/lr_imputation/refs/" + PANEL_NAME + "/RData/ref_package.chr"
vcf_prefix = "results/lr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr"

lr_oneKG_RData, lr_oneKG_vcf_lst, lr_oneKG_vcf_dict = get_vcf_concat_lst(region, ref_prefix, vcf_prefix)

rule lr_imputation_prep_all:
    input:
        bamlist = "results/lr_imputation/bamlist.txt",
        hap = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        legend = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        samples = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

rule lr_imputation_all:
    input:
        RData = [lr_oneKG_RData],
        vcf_regions = [lr_oneKG_vcf_lst],
        vcfs = expand(f"results/lr_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz", chr = chromosome),
        truth = expand("results/lr_imputation/truth/long_read_truth.chr{chr}.vcf.gz", chr = chromosome)

pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['lr_imputation_dir']
lc_vcf_dir = config['lr_lc_vcf_dir']
hc_vcf_dir = config['lr_hc_vcf_dir']

rule lr_imputation_comparison_all:
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