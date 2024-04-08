configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
QUILT_HOME = config["QUILT_HOME"]

rule hla_imputation_preprocess:
    input:
        bam = "data/dedup_bams/{id}.bam"
    output:
        tmp = temp("results/hla/bams/{id}.tmp.bam"),
        chr = "results/hla/bams/{id}.chr6.bam"
    params:
        verbosity = "ERROR",
        sample = "{id}"
    shell: """
        mkdir -p results/hla/bams/

        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {input.bam} \
        -O {output.tmp} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}

        samtools view -o {output.chr} {output.tmp} chr6:25000000-35000000
    """

rule prepare_hla_bamlist:
    input:
        bams = expand("results/hla/bams/{id}.chr6.bam", id = samples_lc)
    output:
        bamlist = "results/hla/imputation/bamlist.txt"
    shell: """
        mkdir -p results/hla/imputation/

        ls results/hla/bams/*.bam > {output.bamlist}
    """

rule hla_imputation:
    input:
        bamlist = rules.prepare_hla_bamlist.output.bamlist,
        ref_dir = "data/hla_ref_panel"
    output:
        vcf = "results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla = tools['quilt_hla']
    shell: """
        mkdir -p results/hla/imputation/genes/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/genes/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.RData} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={QUILT_HOME}hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
    """
