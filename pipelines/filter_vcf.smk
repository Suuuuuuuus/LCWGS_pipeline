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
        af = "data/oneKG_MAFs/oneKG_AF_AFR_chr{chr}.txt"
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

        lcwgsus.save_vcf(lc_af,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )

# The omni5m manifest has col3,4 to be chr and pos
rule prepare_chip_manifest:
    input:
        bpm = config['dense_bpm'],
        egt = config['dense_egt']
    output:
        csv = "data/chip/omni5m/omni5m.csv",
        pos = "data/chip/omni5m/omni5m_sites.tsv"
    resources:
        mem = '30G'
    params:
        picard = tools['picard']
    shell: """
        mkdir -p results/data/chip/omni5m/

        {params.picard} BpmToNormalizationManifestCsv \
        -I {input.bpm} \
        -CF {input.egt} \
        -O {output.csv}

        cut -d ',' -f3,4 {output.csv} | \
        sed 's/,/\t/g' | \
        tail -n +2 > {output.pos}
    """

rule filter_lc_sites:
    input:
        vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{{chr}}.vcf.gz",
        sites = rules.prepare_chip_manifest.output.pos
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af_chip_sites/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        panel = PANEL_NAME,
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        imp_vcf = input.vcf
        chip_sites = input.sites

        lc = lcwgsus.read_vcf(imp_vcf).sort_values(by=['chr', 'pos'])
        metadata = lcwgsus.read_metadata(imp_vcf)

        lc = lc.apply(lcwgsus.convert_to_chip_format, axis = 1)
        
        sites = pd.read_table(chip_sites, sep = '\t', names = common_cols, dtype = {'chr': str, 'pos': int}).drop_duplicates(ignore_index = True)
        sites = sites[sites['chr'] == str(wildcards.chr)]
        sites['chr'] = sites['chr'].astype(int)

        lc_sites = pd.merge(lc, sites, on = common_cols)

        lcwgsus.save_vcf(lc_sites,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af_chip_sites/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )