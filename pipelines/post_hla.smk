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

two_stage_hla_vcf_indir = config["two_stage_hla_vcf_indir"]
two_stage_hla_vcf_outdir = config["two_stage_hla_vcf_outdir"]

def get_two_stage_hla_indir(wildcards):
    ix = two_stage_hla_vcf_outdir.index(wildcards.two_stage_hla_vcf_outdir)
    return two_stage_hla_vcf_indir[ix]

rule filter_chr6_for_hla_imputation:
    input:
        lc_vcf_dir = get_two_stage_hla_indir,
        b38_scaffold = "results/hla/reference/multiEth_sites.b38.vcf.gz"
    output:
        vcf = "{two_stage_hla_vcf_outdir}chr6.vcf.gz",
        tmp_vcf = temp("{two_stage_hla_vcf_outdir}chr6.tmp.vcf.gz"),
        tmp1_vcf = temp("{two_stage_hla_vcf_outdir}chr6.tmp1.vcf.gz"),
        scaffold = temp("{two_stage_hla_vcf_outdir}multiEth_sites.b38.txt")
    params:
        info = config['hla_info_filter']
    shell: """
        mkdir -p {wildcards.two_stage_hla_vcf_outdir}

        bcftools query -f '%CHROM\t%POS\n' {input.b38_scaffold} > {output.scaffold}

        cp -f {input.lc_vcf_dir}/*chr6.*.gz {output.tmp_vcf}
        tabix -f {output.tmp_vcf}

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.tmp1_vcf} {output.tmp_vcf}
        tabix -f {output.tmp1_vcf}

        bcftools view -R {output.scaffold} -m2 -M2 -v snps \
        -Oz -o {output.vcf} {output.tmp1_vcf}
    """