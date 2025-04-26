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
        bam = expand('data/bams/{id}.bam', id = samples_fv)
    output:
        gyp_coverage = 'results/coverage/gyp/gyp_coverage.tsv'
    resources:
        mem = '30G'
    params: coverotron = tools['coverotron']
    localrule: True
    shell: """
        mkdir -p results/coverage/gyp/

        {params.coverotron} -bin 10000 \
        -range chr4:143770000-144240000 \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """