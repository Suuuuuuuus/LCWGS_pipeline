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
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information_expanded.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = '\t')

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
        hlatypes_file = get_hlatypes
    output:
        db_filtered = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_{panel}_only/{gene}.ssv'
    resources:
        mem = '40G'
    threads: 4
    run:
        db = pd.read_csv(input.db, sep = ' ')
        
        if wildcards.gene in HLA_GENES:
            hlatypes = pd.read_csv(input.hlatypes_file, sep = '\t')
            alleles = lcwgsus.extract_unique_two_field_resolution_from_hlatypes(hlatypes, wildcards.gene)

            db_filtered = db[[c for c in db.columns if any(":".join(c.split(':')[:2]) == (f'{wildcards.gene}*{a}') for a in alleles)]]
        else:
            db_filtered = db

        db_filtered.to_csv(output.db_filtered, sep = ' ', index = False, header = True) 

rule save_hla_db_as_fasta:
    input:
        db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_{panel}_only/{gene}.ssv'
    output:
        fasta = temp('/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/{gene}.fasta'),
        fai = '/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/{gene}.index.tsv'
    resources:
        mem = '20G'
    threads: 2
    run:
        os.makedirs(f'/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{wildcards.IPD_IMGT_version}_{wildcards.panel}/', exist_ok=True)
    
        db = pd.read_csv(input.db, sep = ' ')
        df = pd.DataFrame(columns = ['Allele', 'Start', 'End'])
        fake_chr = ''
        num_N = 300

        for c in db.columns:
            refseq = ''.join(db[c].tolist()).replace('.', '').lstrip('*').rstrip('*')
            refseq = refseq.replace('*', 'N')
            df.loc[len(df)] = [c, len(fake_chr), len(fake_chr) + len(refseq) - 1]
            fake_chr = fake_chr + refseq + 'N'*num_N

        lcwgsus.write_db_as_fasta(wildcards.gene, fake_chr, output.fasta)
        df.to_csv(output.fai, sep = '\t', index = False)

rule merge_fastas:
    input:
        fasta = expand('/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/{gene}.fasta', gene = HLA_GENES_ALL_EXPANDED, allow_missing = True)
    output:
        merged_fasta = "/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/HLA.fasta"
    resources:
        mem = '20G'
    threads: 2
    shell: """
        mkdir -p /well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{wildcards.IPD_IMGT_version}_{wildcards.panel}/
        cat {input.fasta} > {output.merged_fasta}

        bwa index {output.merged_fasta}
    """

rule HLA_realignment:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz",
        reference = "/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v3390_oneKG.fasta"
    output:
        bam = "data/realigned_bams/{id}.bam",
        tmp1 = temp("data/realigned_bams/{id}.tmp1.bam")
    resources:
        mem = '40G'
    params: 
        sample = "{id}",
        picard = tools["picard_plus"]
    threads: 8
    shell: """
        mkdir -p data/realigned_bams/

        bwa mem -t {threads} -a {input.reference} {input.fastq1} {input.fastq2} | \
        samtools view -b -o {output.tmp1}
        
        samtools sort -@{threads} -m 1G -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

rule prepare_hla_reference_panel_method_new:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3390.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData", gene = hla_genes)
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

rule hla_imputation_method_new:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_oneKG_only/{gene}.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '40G'
    threads: 4
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.gene}/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.gene}/" \
        --bamlist={output.bamlist} \
        --region={wildcards.gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """