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
'''
rule liftover_multiEth_vcf:
    input:
        vcf = "results/hla/reference/multiEth_sites.b37.vcf.gz",
        reference = "data/references/GRCh38_with_alt.fa",
        chain = "data/ref_panel/b37ToHg38.over.chain",
        dictionary = "data/references/GRCh38_with_alt.dict"
    output:
        tmp_vcf = temp("results/hla/reference/multiEth_sites.b38.tmp.vcf.gz"),
        lifted = "results/hla/reference/multiEth_sites.b38.vcf.gz",
        rejected = "results/hla/reference/multiEth_sites.b38.rejected.vcf.gz"
    resources: mem = '30G'
    threads: 4
    params:
        picard = tools["picard"]
    shell: """
        {params.picard} LiftoverVcf \
        -I {input.vcf} \
        -O {output.tmp_vcf} \
        -CHAIN {input.chain} \
        -REJECT {output.rejected} \
        -WMC true \
        --MAX_RECORDS_IN_RAM 50000 \
        -R {input.reference}

        tabix -f {output.tmp_vcf}

        bcftools view -r chr6 {output.tmp_vcf} | \
        bcftools sort -Oz -o {output.lifted}
        tabix -f {output.lifted}
    """

rule retain_lc_mGenv1_topmed_chr6_chip_sites:
    input:
        lc_vcf = "results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr6.dose.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr6.vcf.gz"
    output:
        filtered_vcf = "results/wip_vcfs/malariaGen_v1_b38/topmed/chip_sites/lc.chr6.vcf.gz",
        site = temp("results/wip_vcfs/malariaGen_v1_b38/topmed/chip_sites/chr6.tsv")
    resources:
        mem = '30G'
    threads: 4
    localrule: True
    shell: """
        mkdir -p results/wip_vcfs/malariaGen_v1_b38/topmed/chip_sites/

        zgrep -v '#' {input.chip_vcf} | cut -f1,2 > {output.site}
        bcftools view -R {output.site} -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """
'''

indir_map2 = dict(zip(config["two_stage_hla_vcf_outdir"], config["two_stage_hla_vcf_indir"]))

def get_two_stage_hla_indir(wildcards):
    return indir_map2.get(wildcards.odir2, "")

rule filter_two_stage_chr6_for_hla_imputation:
    input:
        lc_vcf_dir = get_two_stage_hla_indir,
        b38_scaffold = "results/hla/reference/multiEth_sites.b38.vcf.gz"
    output:
        vcf = "{odir2}chr6.vcf.gz",
        tmp_vcf = temp("{odir2}chr6.tmp.vcf.gz"),
        tmp1_vcf = temp("{odir2}chr6.tmp1.vcf.gz"),
        scaffold = temp("{odir2}multiEth_sites.b38.txt")
    params:
        info = config['hla_info_filter']
    localrule: True
    shell: """
        mkdir -p {wildcards.odir2}

        bcftools query -f '%CHROM\t%POS\n' {input.b38_scaffold} > {output.scaffold}

        cp -f {input.lc_vcf_dir}/*chr6.*.gz {output.tmp_vcf}
        tabix -f {output.tmp_vcf}

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.tmp1_vcf} {output.tmp_vcf}
        tabix -f {output.tmp1_vcf}

        bcftools view -R {output.scaffold} -m2 -M2 -v snps \
        -Oz -o {output.vcf} {output.tmp1_vcf}
    """

indir_map3 = dict(zip(config["three_stage_hla_vcf_outdir"], config["three_stage_hla_vcf_indir"]))

def get_three_stage_hla_indir(wildcards):
    return indir_map3.get(wildcards.odir3, "")

rule filter_three_stage_chr6_for_hla_imputation:
    input:
        lc_vcf_dir = get_three_stage_hla_indir,
        b38_scaffold = "results/hla/reference/multiEth_sites.b38.vcf.gz"
    output:
        vcf = "{odir3}chr6.vcf.gz",
        tmp_vcf = temp("{odir3}chr6.tmp.vcf.gz"),
        scaffold = temp("{odir3}multiEth_sites.b38.txt")
    params:
        info = config['hla_info_filter']
    localrule: True
    shell: """
        mkdir -p {wildcards.odir3}

        bcftools query -f '%CHROM\t%POS\n' {input.b38_scaffold} > {output.scaffold}

        cp -f {input.lc_vcf_dir}/*chr6.*.gz {output.tmp_vcf}
        tabix -f {output.tmp_vcf}

        bcftools view -R {output.scaffold} -m2 -M2 -v snps \
        -Oz -o {output.vcf} {output.tmp_vcf}
    """