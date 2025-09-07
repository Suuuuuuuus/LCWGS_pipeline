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

sv_df_file = 'results/nonahore/denovo/manifest.5K.5percent.tsv'
sv_df = pd.read_csv(sv_df_file, sep = '\t')
denovo_regions = len(sv_df)

rule all:
    input:
        called = expand("results/hla/call/{id}/hla/R1_bestguess_G.txt", id = ['IDT0481'])

rule hla_la_calling:
    input:
        bam = "data/bams/{id}.bam",
        bai = "data/bams/{id}.bam.bai"
    output:
        called = "results/hla/call/{id}/hla/R1_bestguess_G.txt"
    resources:
        mem = '60G'
    threads: 4
    shell: """
        mkdir -p results/hla/call/{wildcards.id}/
        module load Java/17

        HLA-LA.pl \
        --BAM {input.bam} \
        --graph PRG_MHC_GRCh38_withIMGT \
        --workingDir /well/band/users/rbx225/GAMCC/results/hla/call/ \
        --sampleID {wildcards.id} \
        --maxThreads {threads}
    """