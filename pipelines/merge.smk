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

def merge_bam_input(wildcards):
    return expand("data/chunk_bams/tmp/tmp/{id_ary}/{id_ary}.chr{wildcards.chr}.bam", id_ary = test_hc_dict[wildcards.hc])
# def merge_bam_output(wildcards):
#     return expand("data/chunk_bams/tmp/{id}/{id}.chr{chr}.bam", id_ary = test_hc_dict[wildcards.id], chr = chromosome)

# Merging bams
rule merge_bam:
    input:
        bams = merge_bam_input
    output:
        bam = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"), ### Now the wildcards are messed up to avoid intermediate files... Need to come back later
        bai = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai")
    threads: 2
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/chunk_bams/tmp/{wildcards.hc}/
        
        samtools merge {input.bams} | \
        samtools sort -n | \
        samtools fixmate -m | \
        samtools sort | \
        samtools markdup -O BAM -o {output.bam}
        samtools index {output.bam}
    """
# rule index_merge_bam:
#     input:
#         bam = rules.merge_bam.output.bam
#     output:
#         bai = "data/merge_bams/{id}.bam.bai"
#     threads: 8
#     resources:
#         mem = '50G'
#     shell: """
#         samtools index {input.bam}
#     """