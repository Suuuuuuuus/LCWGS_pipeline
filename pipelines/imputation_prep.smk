configfile: "pipelines/config.json"

import os
from os.path import exists
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("scripts")
import lcwgSus
chromosome = [i for i in range(1,23)]

# The followings are global parameters from `activate`:
QUILT_HOME = config["QUILT_HOME"]
ANALYSIS_DIR = config["ANALYSIS_DIR"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]
rmdup=config["rmdup"]

sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)
panels = config["panels"]

rule prepare_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = ids_1x_all)
    output:
        bamlist = "results/imputation/bamlist.txt"
    params:
        threads=1
    shell: """
        mkdir -p {ANALYSIS_DIR}
        if [[ {rmdup} == "True" ]]
        then
            ls data/dedup_bams/*.bam > {output.bamlist}
        else
            ls data/bams/*.bam > {output.bamlist}
        fi
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
        #R -f scripts/make_b38_recomb_map.R --args "./" {RECOMB_POP} {wildcards.chr}
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {ANALYSIS_DIR} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/imputation_refs/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/imputation_refs/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/imputation/refs/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME,
        threads=1
    shell: """
        mkdir -p results/imputation/refs/
        bcftools view --output-file {output.tmp_vcf} --output-type z --min-alleles 2 --max-alleles 2 --types snps {input.vcf}
        tabix {output.tmp_vcf}
        bcftools convert --haplegendsample results/imputation/refs/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/determine_chunks.R"
    output:
        json = "results/imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {ANALYSIS_DIR:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """
