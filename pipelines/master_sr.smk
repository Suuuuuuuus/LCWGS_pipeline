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

chromosome = [i for i in range(1,23)]

read_lengths = ['151-optimal', '151-long', '151-short', '151-real', '300-optimal', '300-long', '300-short', '300-real']
means = ['500', '500', '329', '329', '1000', '1000', '830', '830']
stds = ['13', '13', '13', '13', '50', '50', '200', '200']
error1 = ['0', '0.0024-0.0071', '0', '0.0024-0.0071', '0', '0.0001', '0', '0.0001']
error2 = ['0', '0.0034-0.0105', '0', '0.0034-0.0105', '0', '0.0001', '0', '0.0001']
haplotypes = ['mat', 'pat']
coverage = '0.6'

QUILT_HOME = config["QUILT_HOME"]
sr_analysis_dir = config["sr_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["hc_panel"]

rule short_read_all:
    input:
        fastq1 = expand("data/sr_simulations/{rl}/{rl}.fastq1.gz", rl = read_lengths),
        fastq2 = expand("data/sr_simulations/{rl}/{rl}.fastq2.gz", rl = read_lengths),
        bams = expand("data/sr_bams/{rl}.bam", rl = read_lengths),
        bais = expand("data/sr_bams/{rl}.bam.bai", rl = read_lengths)

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/sr_imputation/regions.json"
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
        file="results/sr_imputation/refs/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file="results/sr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)

rule sr_imputation_prep_all:
    input:
        bamlist = "results/sr_imputation/bamlist.txt",
        # recomb = expand("results/sr_imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        # json = "results/sr_imputation/regions.json",
        # hap = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        # legend = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        # samples = expand(f"results/sr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/sr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/sr_imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule sr_imputation_all:
    input:
        RData = [regions_to_prep],
        vcf_regions = [vcfs_to_impute],
        vcfs = [final_vcfs],
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