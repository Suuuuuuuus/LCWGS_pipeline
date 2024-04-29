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

chromosome = [i for i in range(1,23)]

read_lengths = ['1kb', '2kb', '5kb', '10kb', '20kb']
# read_lengths = ['1kb']
haplotypes = ['mat', 'pat']
# haplotypes = ['mat']
method = 'CCS'
coverage = '0.5'

QUILT_HOME = config["QUILT_HOME"]
lr_analysis_dir = config["lr_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

rule long_read_all:
    input:
        fastqs = expand("data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz", rl = read_lengths, hap = haplotypes),
        bams = expand("data/lr_bams/{rl}.bam", rl = read_lengths),
        bais = expand("data/lr_bams/{rl}.bam.bai", rl = read_lengths)

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/lr_imputation/regions.json"
if os.path.exists(file):
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
        file="results/lr_imputation/refs/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file="results/lr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)

rule lr_imputation_prep_all:
    input:
        bamlist = "results/lr_imputation/bamlist.txt",
        recomb = expand("results/lr_imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        json = "results/lr_imputation/regions.json",
        hap = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        legend = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        samples = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/lr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/lr_imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule lr_imputation_all:
    input:
        RData = [regions_to_prep],
        vcf_regions = [vcfs_to_impute],
        vcfs = [final_vcfs],
        truth = expand("results/lr_imputation/truth/long_read_truth.chr{chr}.vcf.gz", chr = chromosome)

pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['imputation_dir'][-1]
lc_vcf_dir = config['lc_vcf_dir'][-1]
hc_vcf_dir = config['hc_vcf_dir'][-1]

rule lr_imputation_comparison_all:
    input:
        h_report_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_report_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_impacc_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),

        lc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        hc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        afs = expand("{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv", imp_dir = imputation_dir, chr = chromosome),
        
        sumstats = expand("{imp_dir}summary_metrics.tsv", imp_dir = imputation_dir)