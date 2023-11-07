configfile: "pipelines/config.json"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
sys.path.append("scripts")
import lcwgSus

# samples = pd.read_table(config['samples'], header = None, names = ['Code'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)

chromosome = [i for i in range(1,23)]

rule test:
    output:
        imputation_vcf = "results/imputation/tmp/{id}/res.txt"
    params:
        seq_name = lambda wildcards: wildcards.id,
        sample_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Sample_Name'].values,
        chip_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Chip_Name'].values
    shell: """
        echo {params.seq_name}\t{params.sample_name}\t{params.chip_name} > {output.imputation_vcf}
    """
