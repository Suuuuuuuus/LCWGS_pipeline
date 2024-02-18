configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

rule rmdup:
    input:
        bam = "data/bams/{id}.bam"
    output:
        dedup_bam = "data/dedup_bams/{id}.bam"
    resources: mem = '10G'
    shell: """
        samtools rmdup {input.bam} {output.dedup_bam}
    """

rule index_dedup:
    input:
        dedup_bam = rules.rmdup.output.dedup_bam
    output:
        dedup_bai = "data/dedup_bams/{id}.bam.bai"
    resources: mem = '10G'
    shell: """
        samtools index {input.dedup_bam}
    """

rule rmdup_split:
    input:
        bam_chunk = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"
    output:
        tmp1 = temp("data/chunk_bams/{hc}/{hc}.chr{chr}.tmp1.bam"),
        tmp2 = temp("data/chunk_bams/{hc}/{hc}.chr{chr}.tmp2.bam"),
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        dedup_bai_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam.bai",
        metric = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.metrics.txt")
    threads: 2
    resources: mem = '50G'
    shell: """
        samtools sort -o {output.tmp1} {input.bam_chunk}

        picard MarkDuplicates \
        -I {output.tmp1} \
        -O {output.tmp2} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -o {output.dedup_bam_chunk} {output.tmp2}

        samtools index {output.dedup_bam_chunk}
    """
