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
imp_dir = config["imputation_dir"]

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

rule filter_lc_maf:
    input:
        lc_vcf = imp_dir + "vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz",
        af = imp_dir + "vcf/all_samples/af/af.chr{chr}.tsv"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        info = config['info_filter'],
        maf = config['maf_filter'],
        panel = PANEL_NAME
    run:
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        imp_vcf = input.lc_vcf
        af_txt = input.af

        lc = lcwgsus.read_vcf(imp_vcf).sort_values(by=['chr', 'pos'])
        metadata = lcwgsus.read_metadata(imp_vcf)
        af = lcwgsus.read_af(af_txt)

        lc_af = pd.merge(lc, af, on = common_cols)
        lc_af = lc_af[lc_af['MAF'] > float(params.maf)]

        lcwgsus.save_vcf(lc_af,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )
        
        shell("gunzip {output.filtered_vcf}")
        shell("bgzip results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{wildcards.chr}.vcfs")