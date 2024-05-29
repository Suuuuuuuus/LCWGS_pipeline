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

rule prepare_hla_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_lc)
    output:
        bamlist = "results/hla/imputation/bamlist.txt"
    localrule: True
    shell: """
        mkdir -p results/hla/imputation/

        ls data/bams/*.bam > {output.bamlist}
    """

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_ref_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']

rule prepare_hla_reference_panel:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD_IGMT.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand(f"{hla_ref_panel_outdir}HLA{{gene}}fullallelesfilledin.RData", gene = hla_genes)
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt"
    shell: """
        {params.quilt_hla_prep} \
        --outputdir={hla_ref_panel_outdir} \
        --nGen=100 \
        --hla_types_panel={input.hla_types_panel} \
        --ipd_igmt_alignments_zip_file={input.ipd_igmt} \
        --ref_fasta={input.fasta} \
        --refseq_table_file={params.refseq} \
        --full_regionStart=25587319 \
        --full_regionEnd=33629686 \
        --buffer=500000 \
        --region_exclude_file={params.region_exclude_file} \
        --genetic_map_file={input.genetic_map} \
        --reference_haplotype_file={input.hap} \
        --reference_legend_file={input.legend} \
        --reference_sample_file={input.sample} \
        --hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
        --nCores=6
    """

rule hla_imputation:
    input:
        bamlist = "results/hla/imputation/bamlist.txt",
        ref_dir = hla_ref_panel_outdir
    output:
        vcf = "results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict"
    shell: """
        mkdir -p results/hla/imputation/genes/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/genes/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """
'''
rule hla_imputation_alt:
    input:
        bamlist = "results/hla/imputation/bamlist.txt",
        ref_dir = hla_ref_panel_outdir
    output:
        vcf = "results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict" 
    shell: """
        mkdir -p results/hla/imputation/genes/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/genes/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """
'''
