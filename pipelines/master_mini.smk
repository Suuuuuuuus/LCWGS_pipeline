include: "mini.smk"
include: "subsample.smk"

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

samples_lc = read_tsv_as_lst(config['samples_lc'])
panels = config['panels']

concatenate = config['concatenate']

QUILT_HOME = config["QUILT_HOME"]
mini_analysis_dir = config["mini_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

rule mini_all:
    input:
        ss_bams = expand("data/subsampled_bams/{id}_subsampled.bam", id = samples_lc),
        ss_bais = expand("data/subsampled_bams/{id}_subsampled.bam.bai", id = samples_lc),
        bams = expand("data/mini_bams/{id}.bam", id = samples_lc),
        bais = expand("data/mini_bams/{id}.bam.bai", id = samples_lc)
        

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/mini_imputation/regions.json"
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
        file="results/mini_imputation/refs/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file="results/mini_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)

rule mini_imputation_prep_all:
    input:
        bamlist = "results/mini_imputation/bamlist.txt",
        recomb = expand("results/mini_imputation/" + RECOMB_POP + "/" + RECOMB_POP + "-chr{chr}-final.b38.txt.gz", chr = chromosome),
        json = "results/mini_imputation/regions.json",
        hap = expand(f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz", chr = chromosome),
        legend = expand(f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        samples = expand(f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples", chr = chromosome)

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/mini_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/mini_imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

preparations = ['fv', 'mini']

rule mini_imputation_all:
    input:
        RData = [regions_to_prep],
        vcf_regions = [vcfs_to_impute],
        vcfs = [final_vcfs],
        fv = expand(f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/fv/quilt.chr{{chr}}.vcf.gz", chr = chromosome),
        mini = expand(f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/mini/quilt.chr{{chr}}.vcf.gz", chr = chromosome)

rule filter_vcf_all:
    input:
        high_info_vcf = expand(f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info/lc.chr{{chr}}.vcf.gz", chr = chromosome, prep = preparations),
        high_info_high_maf_vcf = expand(f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info_high_af/lc.chr{{chr}}.vcf.gz", chr = chromosome, prep = preparations),
        high_info_high_maf_chip_sites_vcf = expand(f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info_high_af_chip_sites/lc.chr{{chr}}.vcf.gz", chr = chromosome, prep = preparations)

pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['imputation_dir'][-6:]
lc_vcf_dir = config['lc_vcf_dir'][-6:]
hc_vcf_dir = config['hc_vcf_dir'][-6:]

rule mini_imputation_comparison_all:
    input:
        vcfs_all = expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', imp_dir = imputation_dir, chr = chromosome, pair = pair),
        h_report_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_report_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_impacc_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),

        lc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        hc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        afs = expand("{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv", imp_dir = imputation_dir, chr = chromosome),
        
        sumstats = expand("{imp_dir}summary_metrics.tsv", imp_dir = imputation_dir)