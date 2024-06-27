configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

chromosome = [i for i in range(1,23)]

samples_hc = read_tsv_as_lst(config['samples_hc'])

hc_dict = read_tsv_as_dict(samples_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")

rule fix_hc_bam:
    input:
        bam = "data/bams/{id}.bam"
    output:
        tmp1 = temp("data/chunk_bams/tmp/{id}.tmp1.bam")
    threads: 4
    resources: mem = '10G'
    params:
        verbosity = "ERROR",
        sample = "{id}".split("_")[0]
    shell: """
        mkdir -p data/chunk_bams/tmp/

        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {input.bam} \
        -O {output.tmp1} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}
    """

hc_bam_to_concat_dict = {}
for s in samples_hc:
    hc_bam_to_concat_dict[str(s)] = ["data/chunk_bams/tmp/" + k + ".tmp1.bam" for k in hc_dict[str(s)]]

rule merge:
    input:
        bams = lambda wildcards: hc_bam_to_concat_dict[str(wildcards.hc)]
    output:
        bam = "data/merge_bams/{hc}.bam",
        bai = "data/merge_bams/{hc}.bam.bai",
        tmp1 = temp("data/merge_bams/{hc}.tmp1.bam"),
        metric = temp("data/merge_bams/{hc}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/merge_bams/tmp/

        samtools cat -o {output.tmp1} {input.bams}
        samtools sort -@6 -m 1G -T data/merge_bams/tmp/temp{wildcards.hc} -o {output.bam} {output.tmp1}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T data/merge_bams/tmp/temp{wildcards.hc} -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """
