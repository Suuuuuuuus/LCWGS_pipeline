include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
chromosome = [i for i in range(1,23)]

region_file = "data/imputation_accessories/5Mb_chunks_for_coverage.json"
coverotron_output_prefix = 'results/coverage/coverotron/chr'
_, coverotron_lst, coverotron_dict = get_vcf_concat_lst(region_file, '', coverotron_output_prefix, suffix = '.tsv')

rule all:
    input:
        coverage_files = coverotron_lst

rule keep_primary_reads:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/tmp/{id}.bam')
    resources:
        mem = '30G'
    shell: """
        mkdir -p results/coverage/tmp/

        samtools view -h -f 2 -F 2304 {input.bam} |
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """

rule calculate_bin_coverage_per_chunk:
    input:
        bam = expand('results/coverage/tmp/{id}.bam', id = samples_fv)
    output:
        gyp_coverage = 'results/coverage/coverotron/chr{chr}.{regionStart}.{regionEnd}.tsv'
    resources: mem = '30G'
    threads: 4
    params: coverotron = tools['coverotron']
    shell: """
        mkdir -p results/coverage/coverotron/

        {params.coverotron} -bin 1000 \
        -range chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} \
        -threads {threads} \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """
