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

rule lr_alignment:
    input:
        fastq = rules.simulate_reads.output.fastq,
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if config['concatenate'] else config["ref38"]
    output:
        bam = temp("data/lr_bams/tmp/{hap}.{rl}.bam")
    resources:
        mem = '30G'
    threads: 6
    shell: """
        mkdir -p data/lr_bams/tmp/

        pbmm2 align {input.reference} {input.fastq} {output.bam} \
        --rg '@RG\tID:HG02886\tSM:HG02886' \
        --sort -j 4 -J 2
    """

rule lr_clean_bam:
    input:
        bam = rules.lr_alignment.output.bam
    output:
        bam = "data/lr_bams/{hap}.{rl}.bam",
        bai = "data/lr_bams/{hap}.{rl}.bam.bai",
        tmp1 = temp("data/lr_bams/{hap}.{rl}.tmp1.bam"),
        metric = temp("data/lr_bams/{hap}.{rl}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/lr_bams/tmp/{hap}.{rl}/"
    shell: """
        mkdir -p {params.tmpdir}

        samtools index {input.bam}
        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """