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
ANALYSIS_DIR = config["ANALYSIS_DIR"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]
rmdup=config["rmdup"]

samples_lc = read_tsv_as_lst(config['samples_lc'])

rule prepare_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_lc)
    output:
        bamlist = "results/imputation/bamlist.txt"
    shell: """
        mkdir -p {ANALYSIS_DIR}
        if [[ {rmdup} == "True" ]]
        then
            ls data/dedup_bams/*.bam > {output.bamlist}
        else
            ls data/bams/*.bam > {output.bamlist}
        fi
    """
'''
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
        #R -f scripts/make_b38_recomb_map.R --args "./" {RECOMB_POP} {wildcards.chr}
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {ANALYSIS_DIR} {RECOMB_POP} {wildcards.chr}
    """
'''
rule convert_ref:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/imputation/refs/{PANEL_NAME}/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/imputation/refs/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/imputation/refs/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/imputation/refs/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME
    threads: 4
    resources: mem = '30G'
    shell: """
        mkdir -p results/imputation/refs/{params.panel}/

        bcftools norm -m+ {input.vcf} | bcftools view -m2 -M2 -v snps | bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert --haplegendsample results/imputation/refs/{params.panel}/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """
'''
rule determine_chunks:
    input:
        legend = expand(f"results/imputation/refs/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/quilt_accessories/determine_chunks.R"
    output:
        json = "results/imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {ANALYSIS_DIR:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """
'''