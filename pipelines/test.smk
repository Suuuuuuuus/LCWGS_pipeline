configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_test/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_hc = read_tsv_as_lst(config['samples_hc'])
IPD_IMGT_versions = ['3390', '3570']
HLA_GENES_ALL_EXPANDED_PLUS_DRB26789 = HLA_GENES_ALL_EXPANDED + ['DRB26789']

score_diff_in_alignment_genes_ary = [0, 8, 20, 50]
n_mismatches_ary = [1, 3, 5]
weight_ary = ['T', 'F']
extract_dir = "results/hla/imputation/QUILT_HLA_result_method/"

rule all:
    input:
        gyp_coverage = 'results/coverage/gyp/gyp_coverage.tsv'

rule prepare_chip_manifest:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/gyp/tmp/{id}.bam')
    resources:
        mem = '30G'
    shell: """
        mkdir -p results/coverage/gyp/tmp/
        
        samtools view -h -f 2 -F 2304 {input.bam} | 
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """
    
rule calculate_bin_coverage:
    input:
        bam = expand('results/coverage/gyp/tmp/{id}.bam', id = samples_fv[0])
    output:
        gyp_coverage = 'results/coverage/gyp/gyp_coverage.tsv'
    resources:
        mem = '30G'
    params: coverotron = tools['coverotron']
    localrule: True
    shell: """
        {params.coverotron} -bin 10000 \
        -range chr4:143700001-144250000 \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """