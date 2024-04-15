include: "long_read.smk"

#include: "test.smk"
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

chromosome = [i for i in range(1,23)]

# read_lengths = ['1kb', '2kb', '5kb', '10kb', '20kb']
read_lengths = ['1kb', '2kb']
# haplotypes = ['mat', 'pat']
haplotypes = ['mat']
# method = ['clr', 'ccs']
method = 'CCS'
# coverage = '0.5'
coverage = '0.001'

rule long_read_all:
    input:
        fastqs = expand("data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz", rl = read_lengths, hap = haplotypes)