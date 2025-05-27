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

rule all:
    input:
        coverage = 'results/coverage/specific_regions/HBA.tsv'

rule calculate_bin_coverage_per_chunk:
    input:
        bam = expand('data/bams/{id}.bam', id = samples_fv)
    output:
        coverage = 'results/coverage/specific_regions/HBA.tsv'
    localrule: True
    params: 
        coverotron_shared = tools['coverotron_shared'],
        binsize = 300,
        c = 16,
        start = 170000,
        end = 180200
    shell: """
        mkdir -p results/coverage/specific_regions/

        {params.coverotron_shared} -bin {params.binsize} \
        -range chr{params.c}:{params.start}-{params.end} \
        -reads {input.bam} > \
        {output.coverage}
    """