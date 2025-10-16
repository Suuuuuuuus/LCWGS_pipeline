configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

replicates = 100

sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
sv_df = pd.read_csv(sv_df_file, sep = '\t')
eichler_regions = len(sv_df)

sv_df_file = 'results/nonahore/denovo/manifest.5K.5percent.tsv'
sv_df = pd.read_csv(sv_df_file, sep = '\t')
denovo_regions = len(sv_df)

rule all:
    input:
        eichler1 = expand('results/nonahore/eichler/nonahore1/region{eichler}/results.pickle', eichler = [i for i in range(eichler_regions)]),
        eichler2 = expand('results/nonahore/eichler/nonahore2/region{eichler}/results.pickle', eichler = [i for i in range(eichler_regions)]),
        eichler3 = expand('results/nonahore/eichler/nonahore3/region{eichler}/results.pickle', eichler = [i for i in range(eichler_regions)]),
        # denovo = expand('results/nonahore/denovo/region{denovo}/results.pickle', denovo = [i for i in range(denovo_regions)]),
        simulation = expand('results/nonahore/simulate/plausibility/rep{rep}/eval.pickle', rep = [i for i in range(replicates)]),

rule simulate_nonahore:
    output:
        pickle = 'results/nonahore/simulate/plausibility/rep{rep}/eval.pickle'
    threads: 16
    params:
        odir = 'results/nonahore/simulate/plausibility/rep{rep}/'
    script:
        'scripts/simulate_nonahore.py'

rule run_nonahore1_on_eichler:
    input:
        sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
    output:
        pickle = 'results/nonahore/eichler/nonahore1/region{eichler}/results.pickle'
    threads: 4
    params:
        eichler_file = f'{home_dir}recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
        chunk_file = 'data/imputation_accessories/5Mb_chunks_for_coverage.json',
        odir = 'results/nonahore/eichler/nonahore1/region{eichler}/',
        row_ix = '{eichler}'
    script:
        'scripts/run_nonahore1.py'

rule run_nonahore2_on_eichler:
    input:
        sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
    output:
        pickle = 'results/nonahore/eichler/nonahore2/region{eichler}/results.pickle'
    threads: 4
    params:
        eichler_file = f'{home_dir}recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
        chunk_file = 'data/imputation_accessories/5Mb_chunks_for_coverage.json',
        odir = 'results/nonahore/eichler/nonahore2/region{eichler}/',
        row_ix = '{eichler}'
    script:
        'scripts/run_nonahore2.py'

rule run_nonahore3_on_eichler:
    input:
        sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
    output:
        pickle = 'results/nonahore/eichler/nonahore3/region{eichler}/results.pickle'
    threads: 4
    params:
        eichler_file = f'{home_dir}recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
        chunk_file = 'data/imputation_accessories/5Mb_chunks_for_coverage.json',
        odir = 'results/nonahore/eichler/nonahore3/region{eichler}/',
        row_ix = '{eichler}'
    script:
        'scripts/run_nonahore3.py'
        
rule run_nonahore_on_denovo:
    input:
        sv_df_file = 'results/nonahore/denovo/manifest.5K.5percent.tsv'
    output:
        pickle = 'results/nonahore/denovo/region{denovo}/results.pickle'
    threads: 4
    params:
        chunk_file = 'data/imputation_accessories/5Mb_chunks_for_coverage.json',
        odir = 'results/nonahore/denovo/region{denovo}/',
        row_ix = '{denovo}'
    script:
        'scripts/run_nonahore.py'