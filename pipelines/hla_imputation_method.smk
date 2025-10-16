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
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
sys.path.append(f'{homedir}software/QUILT_test/QUILT/Python/')
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

rule filter_hla_db:
    input:
        db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{gene}.ssv',
        hlatypes_file = 'results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv'
    output:
        db_filtered = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_merged_only/{gene}.ssv'
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

rule hla_alignment_matrix:
    input:
        bam = "data/bams/{id}.bam",
        db_files = expand("/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_merged_only/{gene}.ssv", gene = HLA_GENES_ALL_EXPANDED, allow_missing = True)
    output:
        matrix = "results/hla/imputation/WFA_alignments/v{IPD_IMGT_version}/{id}/{gene}/AS_matrix.ssv"
    resources:
        mem = '120G'
    threads: 16
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information_expanded.tsv",
        db_dir = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_merged_only/",
        strict = True
    script:
        "/well/band/users/rbx225/software/QUILT_test/QUILT/Python/hla_calculate_alignment_matrix.py"

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
        matrix = "results/hla/imputation/WFA_alignments/v3390/{id}/{gene}/AS_matrix.ssv"
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '40G'
    threads: 6
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        outputdir = 'results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/',
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p {params.outputdir}
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="{params.outputdir}" \
        --bamlist={output.bamlist} \
        --region={wildcards.gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}quilt.hrc.hla.{wildcards.gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

rule prepare_hla_reference_panel_method_v3570_new:
    input:
        # hla_types_panel = f"{hla_ref_panel_indir}HLA_oneKG_no_ambiguous.txt",
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3570.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/hla{gene}haptypes.RData", gene = hla_genes)
    resources:
        mem = '60G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/"
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

rule hla_imputation_method_v3570_new:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/hla{gene}haptypes.RData",
        matrix = "results/hla/imputation/WFA_alignments/v3570/{id}/{gene}/AS_matrix.ssv"
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '40G'
    threads: 4
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        outputdir = 'results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{gene}/',
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method_v3570/{wildcards.id}/{wildcards.gene}/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="{params.outputdir}" \
        --bamlist={output.bamlist} \
        --region={wildcards.gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

rule prepare_hla_reference_panel_optimal_new:
    input:
        hla_types_panel = f"results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3570.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{gene}haptypes.RData", gene = hla_genes, allow_missing = True),
        exclude_sample_file = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/{id}.tsv")
    resources:
        mem = '50G'
    threads: 5
    conda: "sus1"
    params:
        quilt_hla_prep = tools['quilt_sus_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/"
    shell: """
        echo {wildcards.id} >> {output.exclude_sample_file}

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
        --reference_exclude_samplelist_file={output.exclude_sample_file} \
        --hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
        --nCores=5
    """

rule hla_imputation_optimal_new:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{gene}haptypes.RData",
        matrix = "results/hla/imputation/WFA_alignments/v3570/{id}/{gene}/AS_matrix.ssv"
    output:
        bamfile = temp("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{id}.{gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_optimal/{id}/{gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '20G'
    threads: 2
    params:
        quilt_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        outputdir = "results/hla/imputation/QUILT_HLA_result_optimal/{id}/{gene}/",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/"
    conda: "sus1"
    shell: """
        mkdir -p {params.outputdir}

        echo {input.bam} >> {output.bamfile}
        ulimit -n 50000

        {params.quilt_hla} \
        --outputdir="{params.outputdir}" \
        --bamlist={output.bamfile} \
        --region={wildcards.gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}quilt.hrc.hla.{wildcards.gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

'''

rule hla_alignment_matrix_hc:
    input:
        bam = "data/merge_bams/{hc}.bam",
        db_files = expand("/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_merged_only/{gene}.ssv", gene = HLA_GENES_ALL_EXPANDED, allow_missing = True)
    output:
        matrix = "results/hla/imputation/WFA_alignments/v{IPD_IMGT_version}/{hc}/{gene}/AS_matrix.ssv"
    resources:
        mem = '120G'
    threads: 16
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information_expanded.tsv",
        db_dir = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_merged_only/",
        strict = True
    script:
        "/well/band/users/rbx225/software/QUILT_test/QUILT/Python/hla_calculate_alignment_matrix.py"

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
        bam = "data/bams/{id}.bam",
        reference = "/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/HLA.fasta"
    output:
        bam = "data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}.bam",
        fastq1 = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}_1.fastq"),
        readlist = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}_readlist.txt"),
        tmp1 = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}.tmp1.bam")
    resources:
        mem = '40G'
    params: 
        sample = "{id}",
        picard = tools["picard"]
    threads: 8
    shell: """
        mkdir -p data/realigned_bams/v{wildcards.IPD_IMGT_version}_{wildcards.panel}

        samtools view -h {input.bam} "chr6:25000000-34000000" | \
        samtools sort -n - -o {output.bam}
        
        samtools view {output.bam} | cut -f 1 | uniq > {output.readlist} 

        {params.picard} FilterSamReads \
        -I {input.bam} \
        -O {output.bam} \
        -READ_LIST_FILE {output.readlist} \
        -FILTER includeReadList

        samtools view -F 2052 -h {output.bam} | \
        samtools sort -n - -o {output.tmp1}

        bedtools bamtofastq -i {output.tmp1} -fq {output.fastq1}

        bwa mem -t {threads} -T 0 -h 0 -a {input.reference} {output.fastq1} | \
        samtools view -b -o {output.tmp1}
        
        samtools sort -@{threads} -m 1G -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """
'''