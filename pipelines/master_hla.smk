include: "hla.smk"

include: "auxiliary.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']

rule hla_all:
    input:
        chrs = expand("results/hla/bams/{id}.chr6.bam", id = samples_lc),
        bamlist = "results/hla/imputation/bamlist.txt",
        ref_panel = expand("results/hla/imputation/ref_panel/HLA{gene}fullallelesfilledin.RData", gene = hla_genes),
        hla_imputed = expand("results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes)