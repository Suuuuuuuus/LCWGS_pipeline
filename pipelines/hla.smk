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

rule alignment:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz",
        reference = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        bam = temp("data/bams/tmp/{id}.bam")
    resources:
        mem = '10G'
    threads: 8
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
    """

rule drop_duplicates:
    input:
        bam = rules.alignment.output.bam
    output:
        fixmate = temp("data/bams/tmp/{id}.fixmate.bam"),
        dedup_bam = "data/dedup_bams/{id}.bam"
    resources: mem = '10G'
    shell: """
        samtools view -h {input.bam} chr6 | \
        samtools sort -n - | \
        samtools fixmate -m - - -u | \
        samtools sort - -u | \
        samtools markdup - {output.fixmate}
        samtools index {output.fixmate}

        samtools rmdup {output.fixmate} {output.dedup_bam}
        samtools index {output.dedup_bam}
    """

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

        samtools index {output.tmp}

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

rule prepare_hla_reference_panel:
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
        cd /well/band/users/rbx225/software/QUILT/
        ./QUILT_HLA_prepare_reference.R \
        --outputdir=${reference_package_dir} \
        --nGen=100 \
        --hla_types_panel=${test_dir}/20181129_HLA_types_full_1000_Genomes_Project_panel.txt \
        --ipd_igmt_alignments_zip_file=${test_dir}${ipdigmt_filename} \
        --ref_fasta=${ref_fasta} \
        --refseq_table_file=${refseq_table_file} \
        --full_regionStart=25587319 \
        --full_regionEnd=33629686 \
        --buffer=500000 \
        --region_exclude_file=hla_ancillary_files/hlagenes.txt \
        --genetic_map_file=${inputs_dir}ACB/ACB-chr6-final.b38.txt.gz \
        --reference_haplotype_file=${reference_haplotype_file} \
        --reference_legend_file=${reference_legend_file} \
        --reference_sample_file=${reference_sample_file} \
        --reference_exclude_samples_for_initial_phasing=FALSE \
        --hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
        --nCores=6
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
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={QUILT_HOME}hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
    """
