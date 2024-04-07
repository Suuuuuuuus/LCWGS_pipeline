include: "auxiliary.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_chip = read_tsv_as_lst(config['samples_chip'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
chromosome = [i for i in range(1,23)]

PANEL_NAME = config["PANEL_NAME"]

rule filter_lc_info:
    input:
        lc_vcf = f"results/imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        info = config['info_filter']
    shell: """
        mkdir -p results/wip_vcfs/{PANEL_NAME}/high_info/

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """