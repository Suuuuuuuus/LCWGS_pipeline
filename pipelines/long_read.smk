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

# read_lengths = ['1kb', '2kb', '5kb', '10kb', '20kb']
read_lengths = ['1kb', '2kb']
# haplotypes = ['mat', 'pat']
haplotypes = ['mat']
# method = ['clr', 'ccs']
method = 'CCS'
# coverage = '0.5'
coverage = '0.001'

def get_num_mean_length(wildcards):
    return int(wildcards.rl[:-2])*1000

rule simulate_reads:
    input:
        fasta = "data/lr_fasta/HG02886.{hap}.fa"
    output:
        tmp = temp("data/lr_simulations/{rl}/{hap}.{rl}.fastq"),
        fastq = "data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        model = "data/lr_models/model_qc_" + method.lower(),
        mean_length = get_num_mean_length,
        outdir = "data/lr_simulations/{rl}/",
        prefix = "data/lr_simulations/{rl}/tmp.{hap}.{rl}"
    shell: """
        mkdir -p {params.outdir}

        pbsim --data-type {method} --depth {coverage} --model_qc {params.model} --length-mean {params.mean_length} --prefix {params.prefix} {input.fasta}

        cat {params.prefix}*.fastq > {output.tmp}

        rm {params.prefix}*.fastq
        rm {params.prefix}*.maf
        rm {params.prefix}*.ref

        gzip {output.tmp}
        touch {output.tmp}
    """

rule lr_alignment:
    input:
        fastq = rules.simulate_reads.output.fastq,
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if config['concatenate'] else config["ref38"]
    output:
        bam = temp("data/lr_bams/tmp/{hap}.{rl}.bam")
    resources:
        mem = '30G'
    threads: 6
    shell: """
        mkdir -p data/lr_bams/tmp/

        pbmm2 align {input.reference} {input.fastq} {output.bam} \
        --rg '@RG\tID:HG02886\tSM:HG02886' \
        --sort -j 4 -J 2
    """

rule lr_clean_bam:
    input:
        bam = rules.lr_alignment.output.bam
    output:
        bam = "data/lr_bams/{hap}.{rl}.bam",
        bai = "data/lr_bams/{hap}.{rl}.bam.bai",
        tmp1 = temp("data/lr_bams/{hap}.{rl}.tmp1.bam"),
        metric = temp("data/lr_bams/{hap}.{rl}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/lr_bams/tmp/{hap}.{rl}/"
    shell: """
        mkdir -p {params.tmpdir}

        samtools index {input.bam}
        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

QUILT_HOME = config["QUILT_HOME"]
lr_analysis_dir = config["lr_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

rule prepare_bamlist:
    input:
        bams = expand("data/lr_bams/{hap}.{rl}.bam", rl = read_lengths, hap = haplotypes)
    output:
        bamlist = "results/lr_imputation/bamlist.txt"
    shell: """
        mkdir -p {lr_analysis_dir}
        
        data/lr_bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {lr_analysis_dir} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/lr_imputation/refs/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME,
        threads=1
    shell: """
        mkdir -p results/lr_imputation/refs/
        bcftools view --output-file {output.tmp_vcf} --output-type z --min-alleles 2 --max-alleles 2 --types snps {input.vcf}
        tabix {output.tmp_vcf}
        bcftools convert --haplegendsample results/lr_imputation/refs/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/determine_chunks.R"
    output:
        json = "results/lr_imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {lr_analysis_dir:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/lr_imputation/regions.json"
if os.path.exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file)

rule prepare_ref:
    input:
        json = "results/lr_imputation/regions.json",
        hap = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        recomb = f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/lr_imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem_mb = 30000
    params:
        threads = 8
    shell: """
        mkdir -p results/lr_imputation/refs/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/lr_imputation/refs/RData/other/", \
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
        bamlist = "results/lr_imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/lr_imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
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
        mkdir -p "results/lr_imputation/vcfs/{params.panel}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/lr_imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        bamlist="{input.bamlist}", \
        prepared_reference_filename="{input.RData}", \
        output_filename="{output.vcf}", \
        seed='${{SEED}}')'
    """
    
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

rule concat_quilt_vcf:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/lr_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
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