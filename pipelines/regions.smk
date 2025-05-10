include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
chromosome = [i for i in range(1,23)]
panels = config['panels']

rule all:
    input:
        gyp_coverage = 'results/coverage/gyp/gyp_coverage.tsv'

# regions are take from the start of GYPE - 0.1Mb and GYPA + 0.1Mb rounded to the nearest 10000s
rule prepare_chip_manifest:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/gyp/tmp/{id}.bam')
    resources:
        mem = '30G'
    shell: """
        mkdir -p results/coverage/gyp/tmp/
        
        samtools view -h -f 2 -F 2304 {input.bam} | 
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """
    
rule calculate_bin_coverage:
    input:
        bam = expand('results/coverage/gyp/tmp/{id}.bam', id = samples_fv)
    output:
        gyp_coverage = 'results/coverage/gyp/gyp_coverage.tsv'
    resources:
        mem = '30G'
    params: coverotron = tools['coverotron']
    localrule: True
    shell: """
        {params.coverotron} -bin 10000 \
        -range chr4:143700000-144250000 \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """