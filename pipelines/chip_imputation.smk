configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]
samples_hc = read_tsv_as_lst(config['samples_hc'])

azim_dir = '/well/band/users/rbx225/Azim_chip/'

rule all:
    input:
        gtc_folder = directory(f"{azim_dir}results/gtcs/"),
        tsv = f"{azim_dir}results/chip.tsv",
        vcf = f"{azim_dir}results/chip.vcf.gz",
        lifted = f"{azim_dir}results/chip.b38.vcf.gz",
        rejected = f"{azim_dir}results/chip.b38.rejected.vcf.gz"
    
rule idat2gtc:
    input:
        manifest = f"{azim_dir}manifest/manifest.csv",
        bpm = f"{azim_dir}manifest/manifest.bpm",
        egt = f"{azim_dir}manifest/manifest.egt",
        idats_folder = f"{azim_dir}results/idats/"
    output:
        gtc_folder = directory(f"{azim_dir}results/gtcs/")
    shell: """
        mkdir -p {azim_dir}results/gtcs/
        
        bcftools +idat2gtc \
        --bpm {input.bpm} \
        --egt {input.egt} \
        --idats {input.idats_folder} \
        --output {output.gtc_folder}
    """   

rule gtc2vcf:
    input:
        manifest = f"{azim_dir}manifest/manifest.csv",
        bpm = f"{azim_dir}manifest/manifest.bpm",
        egt = f"{azim_dir}manifest/manifest.egt",
        fasta = "/well/band/users/rbx225/recyclable_files/genomes/b37/human_g1k_v37.fasta",
        gtc_folder = f"{azim_dir}results/gtcs/"
    output:
        tsv = f"{azim_dir}results/chip.tsv",
        vcf = f"{azim_dir}results/chip.vcf.gz"
    resources: mem = '50G'
    shell: """
        mkdir -p {azim_dir}results/tmp/
        
        bcftools +gtc2vcf \
        --do-not-check-bpm \
        --no-version -Ou \
        --csv {input.manifest} \
        --gtcs {input.gtc_folder} \
        --bpm {input.bpm} \
        --egt {input.egt} \
        --fasta-ref {input.fasta} \
        --extra {output.tsv} | \
        bcftools sort -Ou -T {azim_dir}results/tmp/ | \
        bcftools norm --no-version -o {output.vcf} -Oz -c x -f {input.fasta}

        tabix -f {output.vcf}
    """	

rule liftover_chip_vcf:
    input:
        vcf = f"{azim_dir}results/chip.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp1_vcf = temp(f"{azim_dir}results/chip.tmp.b38.vcf.gz"),
        lifted = f"{azim_dir}results/chip.b38.vcf.gz",
        rejected = f"{azim_dir}results/chip.b38.rejected.vcf.gz",
    resources: mem = '80G'
    threads: 4
    params:
        picard = tools["picard_pplus"]
    shell: """
        {params.picard} LiftoverVcf \
        -I {input.vcf} \
        -O {output.tmp1_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        bcftools view {output.tmp1_vcf} | \
        bcftools sort -Oz -o {output.lifted}

        tabix -f {output.lifted}
    """