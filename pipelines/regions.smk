include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
chromosome = [i for i in range(1,23)]

region_file = "data/imputation_accessories/5Mb_chunks_for_coverage.json"
coverotron_output_prefix = 'results/coverage/coverotron/chr'
bam_prefix = 'results/coverage/tmp/chr'
chunk_bams, cov_files_all = get_bam_concat_lst(region_file, samples_fv, bam_prefix, coverotron_output_prefix)

rule all:
    input:
        hba = 'results/coverage/specific_regions/HBA.tsv',
        primary_bams = expand('data/bams/primary/{id}.bam', id = samples_fv),
        # coverage_files = cov_files_all

rule keep_primary_reads:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/tmp/chr{chr}.{regionStart}.{regionEnd}/{id}.bam')
    resources:
        mem = '30G'
    localrule: True
    shell: """
        mkdir -p results/coverage/tmp/chr{wildcards.chr}.{wildcards.regionStart}.{wildcards.regionEnd}/

        samtools view -h -f 2 -F 2304 {input.bam} chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} | \
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """

def get_bams(wildcards):
    return chunk_bams[f'{wildcards.chr}.{wildcards.regionStart}.{wildcards.regionEnd}']

rule calculate_bin_coverage_per_chunk:
    input:
        bam = get_bams
    output:
        gyp_coverage = 'results/coverage/coverotron/chr{chr}.{regionStart}.{regionEnd}.tsv.gz'
    resources: mem = '60G'
    threads: 4
    params: 
        coverotron = tools['coverotron'],
        coverotron_shared = tools['coverotron_shared'],
        bam_dir = "results/coverage/tmp/chr{chr}.{regionStart}.{regionEnd}/"
    shell: """
        mkdir -p results/coverage/coverotron/

        {params.coverotron_shared} -bin 1000 \
        -range chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} \
        -reads {input.bam} > \
        {output.gyp_coverage}

        rm -r {params.bam_dir}*bai
    """

rule keep_primary_reads_specific_region:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = 'data/bams/primary/{id}.bam'
    shell: """
        mkdir -p data/bams/primary/

        samtools view -h -f 2 -F 2304 {input.bam} | \
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """

rule calculate_bin_coverage_specific_region:
    input:
        bam = expand('data/bams/primary/{id}.bam', id = samples_fv)
    output:
        coverage = 'results/coverage/specific_regions/HBA.tsv'
    localrule: True
    params: 
        coverotron_shared = tools['coverotron_shared'],
        binsize = 2000,
        c = 16,
        start = 100000,
        end = 300000
    shell: """
        mkdir -p results/coverage/specific_regions/

        {params.coverotron_shared} -bin {params.binsize} \
        -range chr{params.c}:{params.start}-{params.end} \
        -reads {input.bam} > \
        {output.coverage}
    """