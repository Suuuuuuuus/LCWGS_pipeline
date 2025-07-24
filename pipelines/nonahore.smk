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

sys.path.append('/well/band/users/rbx225/GAMCC/scripts/lcSV/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/Low Coverage Data/gamcc/scripts/lcSV/')
from lcSV import *

replicates = 100

sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
sv_df = pd.read_csv(sv_df_file, sep = '\t')
eichler_regions = len(sv_df)

sv_df_file = 'results/nonahore/denovo/manifest.5K.5percent.tsv'
sv_df = pd.read_csv(sv_df_file, sep = '\t')
denovo_regions = len(sv_df)

rule all:
    input:
        eichler = expand('results/nonahore/eichler/region{eichler}/results.pickle', eichler = [i for i in range(eichler_regions)]),
        # denovo = expand('results/nonahore/denovo/region{denovo}/results.pickle', denovo = [i for i in range(denovo_regions)]),
        simulation = expand('results/nonahore/simulate/plausibility/rep{rep}/eval.pickle', rep = [i for i in range(replicates)]),

rule simulate_nonahore:
    output:
        pickle = 'results/nonahore/simulate/plausibility/rep{rep}/eval.pickle'
    threads: 16
    params:
        odir = 'results/nonahore/simulate/plausibility/rep{rep}/'
    script:
        '/well/band/users/rbx225/GAMCC/scripts/simulate_nonahore.py'

rule run_nonahore_on_eichler:
    input:
        sv_df_file = 'results/nonahore/eichler/manifest.5K.5percent.tsv'
    output:
        pickle = 'results/nonahore/eichler/region{eichler}/results.pickle'
    threads: 4
    params:
        eichler_file = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
        chunk_file = 'data/imputation_accessories/5Mb_chunks_for_coverage.json',
        odir = 'results/nonahore/eichler/region{eichler}/',
        row_ix = '{eichler}'
    script:
        '/well/band/users/rbx225/GAMCC/scripts/run_nonahore.py'

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
        '/well/band/users/rbx225/GAMCC/scripts/run_nonahore.py'