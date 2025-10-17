configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
sys.path.append(f'{home_dir}software/QUILT_test/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
read_lengths = config["sr_read_lengths"]

rule all:
    input:
        uncoverage_rate_all = "results/sr_coverage/uncov_all.tsv",
        uncoverage_rate = expand("results/sr_coverage/uncov/samples/{rl}_uncov.tsv", rl = read_lengths)


rule compute_bedgraph_nozero:
    input:
        bam = "data/sr_bams/{rl}.bam"
    output:
        bedgraph = temp("results/coverage/bedgraphs/{rl}_bedgraph_nozero.bed")
    resources: mem = '50G'
    shell: """
        bedtools genomecov -ibam {input.bam} -bg | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
    """

rule calculate_uncoverage_rate:
    input:
        bam = "data/sr_bams/{rl}.bam",
        bedgraph = "results/coverage/bedgraphs/{rl}_bedgraph_nozero.bed"
    output:
        uncoverage_rate = "results/sr_coverage/uncov/samples/{rl}_uncov.tsv"
    params:
        access_bed = config['access_bed']
    resources:
        mem = '60G'
    threads: 4
    shell: """
        mkdir -p results/sr_coverage/uncov/samples/

        result=$(bedtools coverage -a {params.access_bed} -b {input.bedgraph} -hist | grep all | head -n 1 | cut -f5)
        echo "{wildcards.rl}\t$result" > {output.uncoverage_rate}
    """

rule aggregate_uncoverage_rate:
    input:
        files = expand("results/sr_coverage/uncov/samples/{rl}_uncov.tsv", rl = read_lengths)
    output:
        uncoverage_rate = "results/sr_coverage/uncov_all.tsv"
    localrule: True
    shell: """
        cat {input.files} >> {output.uncoverage_rate}
    """