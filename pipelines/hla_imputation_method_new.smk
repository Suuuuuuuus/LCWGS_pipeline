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
SEED = 42

rule all:
    input:
        db = expand("results/hla/imputation/QUILT_HLA_result_db/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", id = samples_fv, hla_gene = hla_genes),

        imputed_method_v3570 = expand("results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        
        imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),

        extracted_RData = expand("results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/extracted.hla{hla_gene}.RData", hla_gene = hla_genes, id = samples_fv)



rule hla_imputation_db:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData"
    output:
        bamlist = temp("results/hla/imputation/QUILT_HLA_result_db/{id}/{id}.{hla_gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_db/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '20G'
    threads: 2
    conda: "sus"
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/",
        seed = SEED
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_optimal/{wildcards.id}/{wildcards.hla_gene}/
        echo {input.bam} >> {output.bamlist}

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_db/{wildcards.id}/{wildcards.hla_gene}/" \
        --quilt_seed={params.seed} \
        --bamlist={output.bamlist} \
        --region={wildcards.hla_gene} \
        --n_seek_iterations=2 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

# rule hla_imputation_db:
#     input:
#         bamlist = "results/hla/imputation/bamlists_fv/bamlist{num}.txt",
#         ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData"
#     output:
#         imputed = "results/hla/imputation/QUILT_HLA_result_db/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt"
#     resources:
#         mem = '120G'
#     threads: 8
#     conda: "sus"
#     params:
#         quilt_hla = tools['quilt_hla'],
#         fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
#         ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/",
#         seed = SEED
#     shell: """
#         mkdir -p results/hla/imputation/QUILT_HLA_result_db/genes{wildcards.num}/{wildcards.hla_gene}/

#         {params.quilt_hla} \
#         --outputdir="results/hla/imputation/QUILT_HLA_result_db/genes{wildcards.num}/{wildcards.hla_gene}/" \
#         --quilt_seed={params.seed} \
#         --bamlist={input.bamlist} \
#         --region={wildcards.hla_gene} \
#         --prepared_hla_reference_dir={params.ref_dir} \
#         --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
#         --dict_file={params.fa_dict}
#     """

rule hla_imputation_method_v3570_new:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/hla{hla_gene}haptypes.RData",
        matrix = "results/hla/imputation/WFA_alignments/v3570/{id}/{hla_gene}/AS_matrix.ssv"
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{hla_gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '40G'
    threads: 4
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        outputdir = 'results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/',
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/",
        seed = SEED
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method_v3570/{wildcards.id}/{wildcards.hla_gene}/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="{params.outputdir}" \
        --quilt_seed={params.seed} \
        --bamlist={output.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """






rule prepare_hla_reference_panel_optimal_new:
    input:
        hla_types_panel = f"results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv",
        ipd_igmt = f"results/hla/imputation/ref_panel/auxiliary_files/IPD-IMGT-HLA_v3570.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"results/hla/imputation/ref_panel/auxiliary_files/YRI/YRI-chr6-final.b38.txt.gz",
        hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, allow_missing = True),
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
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{hla_gene}haptypes.RData",
        matrix = "results/hla/imputation/WFA_alignments/v3570/{id}/{hla_gene}/AS_matrix.ssv"
    output:
        bamfile = temp("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{id}.{hla_gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '20G'
    threads: 2
    params:
        quilt_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        outputdir = "results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/",
        seed = SEED
    conda: "sus1"
    shell: """
        mkdir -p {params.outputdir}

        echo {input.bam} >> {output.bamfile}
        ulimit -n 50000

        {params.quilt_hla} \
        --outputdir="{params.outputdir}" \
        --bamlist={output.bamfile} \
        --quilt_seed={params.seed} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """


rule extract_QUILT_alignments:
    input:
        RData = "results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/quilt.RDataput.hla{hla_gene}.RData"
    output:
        extracted_RData = "results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/extracted.hla{hla_gene}.RData"
    localrule: True
    script:
        "/well/band/users/rbx225/software/QUILT_test/QUILT/Python/extract_QUILT_alignments.R"