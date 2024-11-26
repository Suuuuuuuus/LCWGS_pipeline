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

rule prepare_hla_bamlist_wip:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_fv)
    output:
        bamlist = expand("results/hla/imputation/bamlists_fv/bamlist{num}.txt", num = bam_numbers),
        bam_all = temp("results/hla/imputation/bamlists_fv/bamlist.txt")
    localrule: True
    params:
        batch = bam_batches
    shell: """
        mkdir -p results/hla/imputation/bamlists_fv/
        for s in {samples_fv}; do
            ls data/bams/$s.bam >> {output.bam_all}
        done

        total_lines=$(wc -l < {output.bam_all})
        lines_per_file=$(( (total_lines + 2) / 3 ))
        split -l $lines_per_file {output.bam_all} results/hla/imputation/prefix_

        i=1
        for file in results/hla/imputation/prefix_*
        do
            mv "$file" "results/hla/imputation/bamlists_fv/bamlist${{i}}.txt"
            i=$((i + 1))
        done
    """

###### Testing db ######

rule prepare_hla_reference_panel_db:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3570.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes)
    resources:
        mem = '60G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/"
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

rule hla_imputation_db:
    input:
        bamlist = "results/hla/imputation/bamlists_fv/bamlist{num}.txt",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData"
    output:
        imputed = "results/hla/imputation/QUILT_HLA_result_db/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '120G'
    threads: 8
    conda: "sus"
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_db/genes{wildcards.num}/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_db/genes{wildcards.num}/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

###### Testing reference panel ######

rule prepare_hla_vcf_merged_ref:
    input:
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr6.vcf.gz",
        tbi = "results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr6.vcf.gz.tbi"
    output:
        tmp_vcf = temp("results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/tmp.chr6.vcf.gz"),
        hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    threads: 6
    resources: mem = '80G'
    params: outdir = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | \
        bcftools view -m2 -M2 -v snps | \
        bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert -h \
        results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6 {output.tmp_vcf}
    """

def convert_idt_to_gam(wildcards):
    ix = samples_fv.index(wildcards.id)
    return samples_fv_gam[ix]

rule prepare_hla_reference_panel_merged_ref:
    input:
        hla_types_panel = f"results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3390.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, allow_missing = True),
        exclude_sample_file = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/{id}.tsv")
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/",
        exclude_sample_gam = convert_idt_to_gam
    shell: """
        echo {params.exclude_sample_gam} >> {output.exclude_sample_file}

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
        --nCores=6
    """

rule hla_imputation_merged_ref:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/HLA{hla_gene}fullallelesfilledin.RData"
    output:
        bamfile = temp("results/hla/imputation/QUILT_HLA_result_merged_ref/{id}/{id}.{hla_gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_merged_ref/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '50G'
    threads: 6
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/"
    conda: "sus"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_merged_ref/{wildcards.id}/{wildcards.hla_gene}/

        echo {input.bam} > {output.bamfile}

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_merged_ref/{wildcards.hla_gene}/" \
        --bamlist={output.bamfile} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

###### Testing imputation method ######

rule prepare_hla_db:
    input:
        ipd_gen_files = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/{hla_gene}_gen.txt'
    output:
        output_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv'
    resources:
        mem = '80G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')

        db = process_db_genfile(wildcards.hla_gene, params.ipd_gen_file_dir, hla_gene_information)
        db.to_csv(output.output_db, sep = ' ', index = False, header = True)

rule prepare_hla_db2:
    input:
        ipd_gen_files = '/well/band/users/rbx225/recyclable_files/hla/alignments/{hla_gene}_gen.txt'
    output:
        output_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/{hla_gene}.ssv'
    resources:
        mem = '80G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')

        db = process_db_genfile(wildcards.hla_gene, params.ipd_gen_file_dir, hla_gene_information)
        db.to_csv(output.output_db, sep = ' ', index = False, header = True)

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

'''
rule prepare_hla_reference_panel_method2:
    input:
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        reads1 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv"),
        reads2 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv"),
        mate_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv"),
        pair_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv")
    resources:
        mem = '120G'
    threads: 8
    params:
        outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}",
        hla_gene_information = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv",
        script = "/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/hla_align.py"
    shell: """
        python {params.script} {wildcards.hla_gene} {input.bam} {params.outdir}
    """
'''

rule hla_imputation_method:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/{hla_gene}.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{hla_gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '60G'
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

###### Optimal run in theory ######

rule prepare_hla_reference_panel_optimal:
    input:
        hla_types_panel = f"results/hla_ref_panel/oneKG_mGenv1/oneKG_mGenv1_HLA_calls.tsv",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v3570.zip",
        fasta = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, allow_missing = True),
        exclude_sample_file = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/{id}.tsv")
    resources:
        mem = '80G'
    threads: 8
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/",
        exclude_sample_gam = convert_idt_to_gam
    shell: """
        echo {params.exclude_sample_gam} >> {output.exclude_sample_file}

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
        --nCores=6
    """

rule hla_imputation_optimal:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/HLA{hla_gene}fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv'
    output:
        bamfile = temp("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{id}.{hla_gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '80G',
        partition = 'long'
    threads: 8
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_optimal/{wildcards.id}/{wildcards.hla_gene}/

        echo {input.bam} > {output.bamfile}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_optimal/{wildcards.hla_gene}/" \
        --bamlist={output.bamfile} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """
