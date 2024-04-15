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

read_lengths = ['1kb', '2kb', '5kb', '10kb', '20kb']
haplotypes = ['mat', 'pat']
# method = ['clr', 'ccs']
method = 'CCS'
# coverage = '0.5'
coverage = '0.001'

def get_num_mean_length(wildcards):
    return int(wildcards.rl[:-2])*1000

rule simulate_reads:
    input:
        fasta = "data/lr_fasta/HG02886.{hap}.fa"
    output:
        tmp = temp("data/lr_simulations/{rl}/{hap}.{rl}.fastq"),
        fastq = "data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        model = "data/lr_models/model_qc_" + method.lower(),
        mean_length = get_num_mean_length,
        outdir = "data/lr_simulations/{rl}/",
        prefix = "data/lr_simulations/{rl}/tmp.{hap}.{rl}"
    shell: """
        mkdir -p {params.outdir}

        pbsim --data-type {method} --depth {coverage} --model_qc {params.model} --length-mean {params.mean_length} --prefix {params.prefix} {input.fasta}

        cat {params.prefix}*.fastq > {output.tmp}

        rm {params.prefix}*.fastq
        rm {params.prefix}*.maf
        rm {params.prefix}*.ref

        gzip {output.tmp}
        touch {output.tmp}
    """