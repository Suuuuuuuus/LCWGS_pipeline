configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

# Adapter trimming
rule trimmomatic_lc:
    input:
        # fastq1 = rules.fastuniq.output.fastq1, 
        # fastq2 = rules.fastuniq.output.fastq2
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        fwd_pair = "data/fastq_cleaned/{id}_1.fastq.gz",
        rev_pair = "data/fastq_cleaned/{id}_2.fastq.gz",
        fwd_unpair = "data/fastq_cleaned/{id}_unpaired_1.fastq.gz",
        rev_unpair = "data/fastq_cleaned/{id}_unpaired_2.fastq.gz"
    params:
        adapters = config['adapter']
    threads: 4
    resources:
        mem = '30G'
    shell: """
        trimmomatic PE \
        {input.fastq1} {input.fastq2} \
        {output.fwd_pair} {output.fwd_unpair} \
        {output.rev_pair} {output.rev_unpair} \
        ILLUMINACLIP:{params.adapters}:2:30:10:2:true MINLEN:50
    """

def sample_hc_fastq1(wildcards):
    sample = (wildcards.id).rsplit('_', 1)[0]
    return "data/fastq/tmp/" + sample + "/" + wildcards.id + "_1.fastq.gz"
def sample_hc_fastq2(wildcards):
    sample = (wildcards.id).rsplit('_', 1)[0]
    return "data/fastq/tmp/" + sample + "/" + wildcards.id + "_2.fastq.gz"

rule trimmomatic_hc:
    input:
        fastq1 = sample_hc_fastq1, 
        fastq2 = sample_hc_fastq2
    output:
        fwd_pair = "data/fastq_cleaned/{id}_1.fastq.gz",
        rev_pair = "data/fastq_cleaned/{id}_2.fastq.gz",
        fwd_unpair = "data/fastq_cleaned/{id}_unpaired_1.fastq.gz",
        rev_unpair = "data/fastq_cleaned/{id}_unpaired_2.fastq.gz"
    params:
        adapters = config['adapter']
    threads: 4
    resources:
        mem = '30G'
    shell: """
        trimmomatic PE \
        {input.fastq1} {input.fastq2} \
        {output.fwd_pair} {output.fwd_unpair} \
        {output.rev_pair} {output.rev_unpair} \
        ILLUMINACLIP:{params.adapters}:2:30:10:2:true MINLEN:50
    """