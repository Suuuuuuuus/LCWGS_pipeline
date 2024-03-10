configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

# samples = pd.read_table(config['samples'], header = None, names = ['Code'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)

chromosome = [str(i) for i in range(1,23)]

ls = ['snps', 'indels']

rule all:
    input:
        expand("tmp/{k}.txt", k = ls)

rule test:
    output:
        vcf = "tmp/{k}.txt"
    params:
    #    mode = "{k}"
        mode = "SNP" if (str("{k}") == "snps") else "INDEL"
    shell: """
        mkdir -p tmp/
        echo {wildcards.k} >> {output.vcf}
        echo {params.mode} >> {output.vcf}
"""
