configfile: "pipelines/config.json"
include: "auxiliary.smk"

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

chromosome = [str(i) for i in range(1,23)]

test_hc = ids_1x_all[:2]
test_hc_dict = read_tsv_as_dict(test_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")

nest = {}
for i in test_hc:
    nest[i] = {}
    for j in chromosome:
        nest[i][j] = ["data/tmp/" + k + "/" + k + ".chr" + str(j) + ".bam" for k in test_hc_dict[i]]


rule test:
    output:
        vcf = "results/tmp/{id}.{chr}.txt"
    shell: """
        echo {wildcards.id}
    """