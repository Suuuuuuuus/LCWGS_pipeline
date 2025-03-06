configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import io
import os
import re
import json
import pandas as pd
import numpy as np
import math
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/')
import lcwgsus

from lcwgsus.variables import *
from hla_phase import *
from hla_align_functions import *
from hla_align import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']

rule prepare_hla_db:
    input:
        ipd_gen_files = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments_v{IPD_IMGT_version}/{gene}_gen.txt'
    output:
        db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{gene}.ssv'
    resources:
        mem = '40G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments_v{IPD_IMGT_version}/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')

        db = process_db_genfile(wildcards.gene, params.ipd_gen_file_dir, hla_gene_information)
        db.to_csv(output.db, sep = ' ', index = False, header = True)

def get_hlatypes(wildcards):
    if wildcards.panel == 'oneKG':
        panel = 'results/hla/imputation/ref_panel/auxiliary_files/20181129_HLA_types_full_1000_Genomes_Project_panel.txt'
    elif wildcards.panel == 'merged':
        panel = 'results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv'
    else:
        panel = ''
    return panel

rule filter_hla_db:
    input:
        db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{gene}.ssv',
        hlatypes = get_hlatypes
    output:
        db_filtered = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_{panel}_only/{gene}.ssv'
    resources:
        mem = '40G'
    threads: 4
    run:
        db = pd.read_csv(input.db, sep = ' ')
        alleles = lcwgsus.extract_unique_two_field_resolution_from_hlatypes(input.hlatypes, wildcards.gene)

        db_filtered = db[[c for c in db.columns if any(":".join(c.split(':')[:2]) == (f'{wildcards.gene}*{a}') for a in alleles)]]
        db_filtered.to_csv(output.db_filtered, sep = ' ', index = False, header = True)        

rule prepare_hla_reference_panel_method:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3390.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes)
    resources:
        mem = '60G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    shell: """
        {params.quilt_hla_prep} \
        --outputdir={params.hla_ref_panel_outdir} \
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

rule hla_imputation_method:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/{hla_gene}.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{hla_gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '40G'
    threads: 4
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.hla_gene}/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.hla_gene}/" \
        --bamlist={output.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """