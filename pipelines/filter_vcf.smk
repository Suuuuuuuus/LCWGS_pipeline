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
        lc_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info/lc.chr{{chr}}.vcf.gz",
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        info = config['info_filter'],
        maf = config['maf_filter'],
        panel = PANEL_NAME,
        chrom = "{chr}"
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
        lc_af = lc_af.drop(columns = 'MAF')
        lc_af = lc_af.apply(lcwgsus.convert_to_chip_format, axis = 1)

        lcwgsus.save_vcf(lc_af,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )
        
        shell("""
            gunzip results/wip_vcfs/{panel}/high_info_high_af/lc.chr{c}.vcf.gz; bgzip results/wip_vcfs/{panel}/high_info_high_af/lc.chr{c}.vcf
        """.format(panel = params.panel, c = params.chrom))