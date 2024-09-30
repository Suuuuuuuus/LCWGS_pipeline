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
samples_fv = read_tsv_as_lst("data/sample_tsvs/fv_idt_names.tsv")
chromosome = [i for i in range(1,23)]
QUILT_HOME = config["QUILT_HOME"]
NGEN = config["NGEN"]

rule all:
    input:
        hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples",
        RData = "results/hla_tests/quilt.hrc.hla.all.haplotypes.RData",
        bam_all = "results/hla_tests/bamlist.txt"

rule subset_vcf_to_chr6:
    input:
        vcf = "results/wip_vcfs/oneKG/vanilla/high_info_high_af_high_conf/lc.chr6.vcf.gz"
    output:
        tmp_vcf = temp("results/hla_tests/gamcc_vcf/fv.chr6.vcf.gz"),
        hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples"
    threads: 4
    resources: mem = '30G'
    params: 
        outdir = "results/hla_tests/gamcc_vcf/",
        fv = "data/sample_tsvs/fv_gm_names.tsv"
    shell: """
        mkdir -p {params.outdir}

        bcftools view -S {params.fv} {input.vcf}| \
        bcftools norm -m+ | \
        bcftools view -m2 -M2 -v snps | \
        
        bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert -h \
        {params.outdir}fv.chr6 {output.tmp_vcf}

        sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' {output.samples}
    """


rule prepare_hla_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_fv)
    output:
        bam_all = "results/hla_tests/bamlist.txt"
    localrule: True
    shell: """
        mkdir -p results/hla_tests/
        ls {input.bams} > {output.bam_all}
    """

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']

rule prepare_ref:
    input:
        hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples",
        genetic_map = "data/imputation_accessories/maps/YRI-chr6-final.b38.txt"
    output:
        RData = "results/hla_tests/quilt.hrc.hla.all.haplotypes.RData"
    resources:
        mem = '30G'
    threads: 4
    params:
        outputdir = "results/hla_tests/prepared_ref/"
    shell: """
        mkdir -p {params.outputdir}
        R -e 'library("data.table"); library("QUILT");
        QUILT_prepare_reference(
        outputdir = "{params.outputdir}",
        nGen = {NGEN},
        chr = "chr6",
        regionStart = 25587319,
        regionEnd = 33629686,
        buffer = 500000,
        reference_haplotype_file = "{input.hap}",
        reference_legend_file = "{input.legend}",
        reference_sample_file = "{input.samples}",
        genetic_map_file = "{input.genetic_map}",
        reference_exclude_samplelist_file = "",
        output_file = "{output.RData}")'
    """

'''
rule hla_imputation:
    input:
        bamlist = "results/hla/imputation/bamlists/bamlist{num}.txt",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_v{IPD_IMGT_version}/"
    output:
        imputed = "results/hla/imputation/QUILT_HLA_result_v{IPD_IMGT_version}/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '50G'
    threads: 6
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_v{IPD_IMGT_version}/genes{wildcards.num}/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_v{wildcards.IPD_IMGT_version}/genes{wildcards.num}/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """
'''
