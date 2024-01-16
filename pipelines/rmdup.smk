configfile: "pipelines/config.json"

from os.path import exists
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
    resources: mem = '100G'
    shell: """
        samtools rmdup {input.bam} {output.dedup_bam}
    """

rule index_dedup:
    input:
        dedup_bam = rules.rmdup.output.dedup_bam
    output:
        dedup_bai = "data/dedup_bams/{id}.bam.bai"
    resources: mem = '100G'
    shell: """
        samtools index {input.dedup_bam}
    """
