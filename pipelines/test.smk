configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

imp_dir = config['imputation_dir']

rule all:
    input:
        f"{imp_dir}vcf/all_samples/lc_vcf/test.txt"

rule test:
    output:
        f"{imp_dir}vcf/all_samples/lc_vcf/test.txt"
    shell: """
        mkdir -p {imp_dir}vcf/
        mkdir -p {imp_dir}impacc/
        mkdir -p {imp_dir}graphs/
        mkdir -p {imp_dir}vcf/all_samples/lc_vcf/
        mkdir -p {imp_dir}vcf/all_samples/hc_vcf/
        mkdir -p {imp_dir}vcf/by_cc/lc_vcf/
        mkdir -p {imp_dir}vcf/by_cc/hc_vcf/
        mkdir -p {imp_dir}vcf/by_eth/lc_vcf/
        mkdir -p {imp_dir}vcf/by_eth/hc_vcf/

        touch {imp_dir}vcf/all_samples/lc_vcf/test.txt
    """