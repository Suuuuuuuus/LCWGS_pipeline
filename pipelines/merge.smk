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

samples_hc = read_tsv_as_lst(config['samples_hc'])

rule fix_bam:
    input:
        bam = "data/bams/{id}.bam"
    output:
        tmp1 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.tmp1.bam"),
        tmp2 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.tmp2.bam")
    threads: 4
    resources: mem = '10G'
    params:
        verbosity = "ERROR",
        sample = "{id}".split("_")[0]
    shell: """
        mkdir -p data/chunk_bams/tmp/tmp/{wildcards.id}/
        samtools sort -o {output.tmp1} {input.bam}

        samtools index {output.tmp1}

        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {output.tmp1} \
        -O {output.tmp2} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}

        picard FixMateInformation -I {output.tmp2}
    """

hc_dict = read_tsv_as_dict(samples_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")

nest = {}
for i in samples_hc:
    nest[str(i)] = ["data/chunk_bams/tmp/tmp/" + k + "/" + k + ".tmp2.bam" for k in hc_dict[str(i)]]

rule merge:
    input:
        bams = lambda wildcards: nest[str(wildcards.hc)]
    output:
        bam = "data/merge_bams/tmp/{hc}.bam",
        bai = "data/merge_bams/tmp/{hc}.bam.bai",
        tmp1 = temp("data/merge_bams/tmp/{hc}.tmp1.bam"),
        metric = temp("data/merge_bams/tmp/{hc}.metrics.txt")
    threads: 2
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/chunk_bams/tmp/{wildcards.hc}/
        mkdir -p data/merge_bams/tmp/
        samtools cat -o {output.tmp1} {input.bams}
        samtools sort -o {output.bam} {output.tmp1}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -o {output.bam} {output.tmp1}
        
        samtools index {output.bam}
    """