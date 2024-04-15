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
        mean_length = int("{rl}"[:-2])*1000,
        outdir = "data/lr_simulations/{rl}/"
    shell: """
        mkdir -p {params.outdir}

        pbsim --data-type {method} --depth {coverage} --model_qc {params.model} {input.fasta} --length-mean {params.mean_length} --prefix {params.outdir}tmp.{wildcards.hap}.{wildcards.rl}

        cat {params.outdir}tmp.{wildcards.hap}.{wildcards.rl}*.fastq > {output.tmp}

        rm {params.outdir}tmp.{wildcards.hap}.{wildcards.rl}*.fastq
        rm {params.outdir}tmp.{wildcards.hap}.{wildcards.rl}*.maf
        rm {params.outdir}tmp.{wildcards.hap}.{wildcards.rl}*.ref

        gzip {output.tmp}
    """