configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

chromosome = [i for i in range(1,23)]

sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
test_hc = ids_1x_all[:2]
test_hc_dict = read_tsv_as_dict(test_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")

nest = {}
for i in test_hc:
    nest[i] = {}
    for j in chromosome:
        nest[str(i)][str(j)] = ["data/chunk_bams/tmp/tmp/" + k + "/" + k + ".chr" + str(j) + ".bam" for k in test_hc_dict[str(i)]]

''' For testing not to split bams into chromosomes

# Merging bams
rule merge_splited_bam:
    input:
        bams = lambda wildcards: nest[str(wildcards.hc)][str(wildcards.chr)]
    output:
        # bam = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"), ### Now the wildcards are messed up to avoid intermediate files... Need to come back later
        # bai = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai")
        bam = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"
#        bai = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai"
    threads: 2
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/chunk_bams/tmp/{wildcards.hc}/
        samtools cat -o {output.bam} {input.bams}
    """

'''