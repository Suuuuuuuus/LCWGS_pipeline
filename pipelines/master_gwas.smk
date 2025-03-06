include: "gwas.smk"
include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

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
import lcwgsus
from lcwgsus.variables import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gm = read_tsv_as_lst('data/sample_tsvs/fv_gm_names.tsv')
samples_oneKG = read_tsv_as_lst("/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/samples_to_phase.tsv")
chromosome = [i for i in range(1,23)]

# models = ['mlr', 'lr1', 'lr2', 'lr3']
models = ['lr1']
sources = ['lc', 'chip']

rule gwas_all:
    input:
        concat_chip_sites_lc = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz",
        PCs_lc = "results/wip_vcfs/oneKG/vanilla/chip_sites/PCs.eigenvec",
        blood_groups_vcf = expand("results/gwas/mGenv1_topmed/{model}/{source}/vcf/blood_groups.vcf.gz", model = models, source = sources),
        gamcc_gen_samples = expand("results/gwas/mGenv1_topmed/{model}/{source}/vcf/blood_groups.sample", model = models, source = sources),
        gamcc_gen = expand("results/gwas/mGenv1_topmed/{model}/{source}/vcf/blood_groups.gen", model = models, source = sources),
        result = expand("results/gwas/mGenv1_topmed/{model}/{source}/results/stats.out", model = models, source = sources),

        PCs_chip = "results/chip/PCs/PCs.eigenvec",
