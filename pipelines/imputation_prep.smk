configfile: "pipelines/config.json"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
chromosome = [i for i in range(1,23)]

QUILT_HOME = config["QUILT_HOME"]
RECOMB_POP = config["RECOMB_POP"]
NGEN = config["NGEN"]
WINDOWSIZE = config["WINDOWSIZE"]
BUFFER = config["BUFFER"]
panels = config['panels']

samples_lc = read_tsv_as_lst(config['samples_lc'])

rule prepare_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_lc)
    output:
        bamlist = "results/imputation/bamlist.txt"
    shell: """
        mkdir -p results/imputation/

        ls data/bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args results/imputation/ {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = "data/ref_panel/{panel}/{panel}.chr{chr}.vcf.gz",
        tbi = "data/ref_panel/{panel}/{panel}.chr{chr}.vcf.gz.tbi"
    output:
        tmp_vcf = temp("results/imputation/refs/{panel}/tmp.{panel}.chr{chr}.vcf.gz"),
        hap = "results/imputation/refs/{panel}/{panel}.chr{chr}.hap.gz",
        legend = "results/imputation/refs/{panel}/{panel}.chr{chr}.legend.gz",
        samples = "results/imputation/refs/{panel}/{panel}.chr{chr}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    threads: 4
    resources: mem = '30G'
    params: outdir = "results/imputation/refs/{panel}/"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | bcftools view -m2 -M2 -v snps | bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert -h \
        {params.outdir}{wildcards.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand("results/imputation/refs/{panel}/{panel}.chr{chr}.legend.gz", chr = chromosome, allow_missing = True)
    output:
        json = "results/imputation/refs/{panel}/regions.json"
    resources: mem = '10G'
    params:
        analysis_dir = "results/imputation/",
        code = "scripts/quilt_accessories/determine_chunks.R"
    shell: """
        Rscript {params.code} {params.analysis_dir} {WINDOWSIZE} {BUFFER} {wildcards.panel}
    """