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

QUILT_HOME = config["QUILT_HOME"]
mini_analysis_dir = config["mini_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

samples_lc = read_tsv_as_lst(config['samples_lc'])

rule prepare_bamlist:
    input:
        bams = expand("data/subsampled_bams/{id}_subsampled.bam", id = samples_lc)
    output:
        bamlist = "results/mini_imputation/bamlist.txt"
    shell: """
        mkdir -p {mini_analysis_dir}
        
        ls data/subsampled_bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {mini_analysis_dir} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/mini_imputation/refs/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME,
        threads=1
    shell: """
        mkdir -p results/mini_imputation/refs/
        bcftools norm -m+ {input.vcf} | bcftools view -Oz -o {output.tmp_vcf} -m2 -M2 -v snps

        tabix {output.tmp_vcf}

        bcftools convert --haplegendsample results/mini_imputation/refs/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/determine_chunks.R"
    output:
        json = "results/mini_imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {mini_analysis_dir:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """

rule prepare_ref:
    input:
        json = "results/mini_imputation/regions.json",
        hap = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        recomb = f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/mini_imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem_mb = 30000
    params:
        threads = 8
    shell: """
        mkdir -p results/mini_imputation/refs/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/mini_imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        nGen={NGEN}, \
        reference_haplotype_file="{input.hap}" ,\
        reference_legend_file="{input.legend}", \
        genetic_map_file="{input.recomb}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        output_file="{output.RData}")'
    """

rule quilt_imputation:
    input:
        bamlist = "results/mini_imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/mini_imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem_mb = 30000
    resources:
        mem_mb = 30000
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    params:
        panel = PANEL_NAME
    shell: """
        ## set a seed here, randomly, so can try to reproduce if it fails
        SEED=`echo $RANDOM`
        mkdir -p "results/mini_imputation/vcfs/{params.panel}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/mini_imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        bamlist="{input.bamlist}", \
        prepared_reference_filename="{input.RData}", \
        output_filename="{output.vcf}", \
        seed='${{SEED}}')'
    """

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/mini_imputation/regions.json"
if os.path.exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file)

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

rule concat_quilt_vcf:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/mini_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
    resources:
        mem_mb = 30000
    params:
        threads = 1,
        input_string=get_input_vcfs_as_string
        # rename_samples = config["rename_samples"],
        # rename_samples_file = config["rename_samples_file"]
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        if [ -e {output.vcf}.temp1.vcf.gz ]; then
            rm {output.vcf}.temp1.vcf.gz
        fi
        if [ -e {output.vcf}.temp2.vcf.gz ]; then
            rm {output.vcf}.temp2.vcf.gz
        fi

        bcftools concat \
        --ligate-force \
        --output-type z \
        --output {output.vcf}.temp1.vcf.gz \
        {params.input_string}

        gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS\t[%GT:%GP:%DS\t]\n' {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf
        bcftools sort -Oz -o {output.vcf} {output.vcf}.temp2.vcf
        tabix {output.vcf}
        rm {output.vcf}.temp*
    """
