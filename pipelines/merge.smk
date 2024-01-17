configfile: "pipelines/config.json"
include: "auxiliary.smk"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
test_hc = ids_1x_all[:2]
test_hc_dict = read_tsv_as_dict(test_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")

def sample_hc_to_sample_hc_lst(wildcards):
    return expand("data/bams/{id_ary}.bam", id_ary = test_hc_dict[wildcards.id])

# Merging bams
rule merge_bam:
    input:
        bams = sample_hc_to_sample_hc_lst
    output:
        bam = "data/merge_bams/{id}.bam"
    threads: 8
    resources:
        mem = '50G'
    shell: """
        samtools merge -O BAM {output.bam} {input.bams}
    """
rule index_merge_bam:
    input:
        bam = rules.merge_bam.output.bam
    output:
        bai = "data/merge_bams/{id}.bam.bai"
    threads: 8
    resources:
        mem = '50G'
    shell: """
        samtools index {input.bam}
    """