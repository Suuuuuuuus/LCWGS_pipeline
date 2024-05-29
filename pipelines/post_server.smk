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

two_stage_vcf_indir = config["two_stage_vcf_indir"]
two_stage_vcf_outdir = config["two_stage_vcf_outdir"]

def get_indir_vcf(wildcards):
    ix = two_stage_vcf_outdir.index(wildcards.two_stage_vcf_outdir)
    return two_stage_vcf_indir[ix] + "chr6.dose.vcf.gz"

rule further_filter_vcf:
    input:
        vcf = get_indir_vcf
    output:
        vcf = temp("{two_stage_vcf_outdir}chr6.tmp.vcf.gz")
    resources:
        mem = '30G'
    threads: 4
    shell:"""
        mkdir -p {wildcards.two_stage_vcf_outdir}

        bcftools view -r chr6:25000000-34000000 -Oz -o {output.vcf} {input.vcf}
    """

rule convert_chr6_to_chip_form:
    input:
        vcf = "{two_stage_vcf_outdir}chr6.tmp.vcf.gz"
    output:
        vcf = "{two_stage_vcf_outdir}chr6.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    run:
        imp_vcf = input.vcf
        lc = lcwgsus.read_vcf(imp_vcf)
        metadata = lcwgsus.read_metadata(imp_vcf)

        lc = lc.apply(lcwgsus.convert_to_chip_format, axis = 1)

        lcwgsus.save_vcf(lc,
             metadata,
             prefix='chr',
             outdir=wildcards.two_stage_vcf_outdir,
             save_name="chr6.vcf.gz"
             )
        lcwgsus.rezip_vcf(output.vcf)