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
sys.path.append('/well/band/users/rbx225/software/QUILT_test/QUILT/Python/')
import lcwgsus

from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']
HLA_GENES_ALL_EXPANDED_PLUS_DRB26789 = HLA_GENES_ALL_EXPANDED + ['DRB26789']

rule all:
    input:
        db = expand('results/hla/align/dbs/v{IPD_IMGT_version}/{gene}.ssv', gene = HLA_GENES_ALL_EXPANDED, IPD_IMGT_version = IPD_IMGT_versions[0]),
        fasta = expand('results/hla/align/fastas_MSA/v{IPD_IMGT_version}/{gene}/{gene}.fasta', gene = HLA_GENES_ALL_EXPANDED, IPD_IMGT_version = IPD_IMGT_versions[0]),
        fasta_DRB26789 = expand('results/hla/align/fastas_MSA/v{IPD_IMGT_version}/DRB26789/DRB26789.fasta', IPD_IMGT_version = IPD_IMGT_versions[0]),
        
        graph_gfa = expand('results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}.gfa', gene = HLA_GENES_ALL_EXPANDED_PLUS_DRB26789, IPD_IMGT_version = IPD_IMGT_versions[0]),
        # index = expand('results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}.giraffe.gbz', gene = HLA_GENES_ALL_EXPANDED_PLUS_DRB26789, IPD_IMGT_version = IPD_IMGT_versions[0]),

rule prepare_hla_db_new:
    input:
        ipd_gen_files = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments_v{IPD_IMGT_version}/{gene}_gen.txt'
    output:
        db = 'results/hla/align/dbs/v{IPD_IMGT_version}/{gene}.ssv'
    resources:
        mem = '40G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments_v{IPD_IMGT_version}/'
    run:
        hla_gene_information = pd.read_csv(HLA_GENE_INFORMATION_FILE_EXPANDED, sep = ' ')

        db = process_db_genfile(wildcards.gene, params.ipd_gen_file_dir, hla_gene_information)
        db.to_csv(output.db, sep = ' ', index = False, header = True)

rule save_hla_db_as_fasta_per_allele:
    input:
        db = 'results/hla/align/dbs/v{IPD_IMGT_version}/{gene}.ssv'
    output:
        fasta = 'results/hla/align/fastas_MSA/v{IPD_IMGT_version}/{gene}/{gene}.fasta'
    resources:
        mem = '20G'
    threads: 2
    run:
        os.makedirs(f'results/hla/align/fastas_MSA/v{wildcards.IPD_IMGT_version}/{wildcards.gene}', exist_ok=True)

        db = pd.read_csv(input.db, sep = ' ')
        lcwgsus.write_db_as_fasta_per_allele(db, output.fasta)

rule save_DRB26789:
    input:
        fasta = '/well/band/users/rbx225/recyclable_files/hla_reference_files/DRB26789_MSA.fasta'
    output:
        fasta = 'results/hla/align/fastas_MSA/v{IPD_IMGT_version}/DRB26789/DRB26789.fasta'
    resources:
        mem = '20G'
    threads: 2
    shell: """
        cp {input.fasta} {output.fasta}
    """

rule construct_graph:
    input: 
        fasta = 'results/hla/align/fastas_MSA/v{IPD_IMGT_version}/{gene}/{gene}.fasta'
    output: 
        fasta = temp('results/hla/align/fastas_MSA/v{IPD_IMGT_version}/{gene}/{gene}.tmp.fasta'),
        graph_vg = temp('results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}.vg'),
        graph_gfa = 'results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}.gfa',
        gbz_index = 'results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}.gbz',

    params:
        prefix = 'results/hla/align/graphs/v{IPD_IMGT_version}/{gene}/{gene}'
    resources:
        mem = '60G'
    threads: 6
    shell: """
        mkdir -p results/hla/align/graphs/v{wildcards.IPD_IMGT_version}/{wildcards.gene}/

        sed '/^>/!s/-/N/g' {input.fasta} > {output.fasta}
        bwa index {output.fasta}

        vg construct -r {output.fasta} > {output.graph_vg}
        vg view {output.graph_vg} > {output.graph_gfa}

        vg autoindex --workflow giraffe -g {output.graph_gfa} -p {params.prefix}
    """
'''
rule align_reads:
    input:
        fq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        gg="hla_graph.gg"
    output: "{sample}.gam"
    shell:
        "vg giraffe -Z {input.gg} -f {input.fastq} -o GAM > {output}"

rule convert_to_bam:
    input: "{sample}.gam"
    output: "{sample}_mapped.bam"
    shell:
        "vg view -a {input} -X | samtools view -bS - | samtools sort -o {output} && "
        "samtools index {output}"
'''