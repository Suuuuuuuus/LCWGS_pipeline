configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append("/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/")
import lcwgsus
from lcwgsus.variables import *
from hla_phase import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_fv = read_tsv_as_lst("data/sample_tsvs/fv_idt_names.tsv")
chromosome = [i for i in range(1,23)]
QUILT_HOME = config["QUILT_HOME"]
NGEN = config["NGEN"]
RECOMB_POP = config["RECOMB_POP"]
studies = ['1KG', 'GAMCC']
to_merge = ['gamcc', 'oneKG']

rule pre_prepare_merge_GAMCC_vcf:
    input:
        vcf = "results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.dose.vcf.gz"
    output:
        fv_vcf = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf.gz",
        tmp_vcf = temp("results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf")
    resources: mem = '30G'
    threads: 4
    params: 
        outdir = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/",
        fv_gm_names = "data/sample_tsvs/fv_gm_names.tsv",
        gm_to_gam = "data/rename_tsvs/fv_gm_to_gam.ssv"
    shell: """
        mkdir -p {params.outdir}

        bcftools view -S {params.fv_gm_names} {input.vcf} | \
        bcftools reheader -s {params.gm_to_gam} -o {output.tmp_vcf}

        bgzip {output.tmp_vcf}
        touch {output.tmp_vcf}
    """

def get_vcf_from_to_merge(wildcards):
    if wildcards.to_merge == 'oneKG':
        return "data/ref_panel/oneKG_30x/oneKG.chr" + wildcards.chr + ".vcf.gz"
    elif wildcards.to_merge == 'gamcc':
        return "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr" + wildcards.chr + ".vcf.gz"
    else:
        return ""

rule prepare_merge_1KG_GAMCC_vcf:
    input:
        vcf = get_vcf_from_to_merge
    output:
        haps = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.hap"),
        legend = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.legend")
    resources: mem = '30G'
    threads: 4
    params: outdir = "results/hla_ref_panel/oneKG_mGenv1/tmp/"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | \
        bcftools view -m2 -M2 -v snps | \
        bcftools sort | \
        bcftools convert -h {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}

        gunzip -f {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.hap.gz
        gunzip -f {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.legend.gz
    """

rule prepare_merge_1KG_GAMCC_sample:
    input:
        mg = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr6.vcf.gz",
        oneKG = "data/ref_panel/oneKG_30x/oneKG.chr6.vcf.gz"
    output:
        tmp_sample = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/tmp.sample"),
        sample = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/merged.samples")
    resources: mem = '30G'
    threads: 1
    params: outdir = "results/hla_ref_panel/oneKG_mGenv1/tmp/"
    shell: """
        mkdir -p {params.outdir}

        echo "sample population group sex" >> {output.sample}

        bcftools query -l {input.oneKG} >> {output.tmp_sample}
        bcftools query -l {input.mg} >> {output.tmp_sample}

        for i in $(cat {output.tmp_sample}); do
            echo $i $i $i 2 >> {output.sample}
        done
    """

rule merge_1KG_GAMCC_per_chunk:
    input:
        haps = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.hap", to_merge = to_merge, allow_missing = True),
        legends = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.legend", to_merge = to_merge, allow_missing = True),
        gen_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr{{chr}}-final.b38.txt",
        sample = rules.prepare_merge_1KG_GAMCC_sample.output.sample
    output:
        haps = temp("results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.hap"),
        legend = temp("results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.legend"),
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.vcf.gz",
    resources: mem = '100G'
    threads: 8
    params: 
        impute2 = tools['impute2'],
        outdir = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/",
        output_prefix = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}",
        mGen_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr{chr}.hap",
        mGen_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr{chr}.legend",
        oneKG_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr{chr}.hap",
        oneKG_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr{chr}.legend",
    shell: """
        mkdir -p {params.outdir}

       {params.impute2} \
        -merge_ref_panels_output_ref {params.output_prefix}.tmp \
        -m {input.gen_map} \
        -h {params.oneKG_haps} \
           {params.mGen_haps} \
        -l {params.oneKG_legend} \
           {params.mGen_legend} \
        -int {wildcards.regionStart} {wildcards.regionEnd} \
        -k_hap 2018 420 \
        -Ne 20000

        awk -F ' ' 'NR==1 {{print; next}} {{$1 = "chr{wildcards.chr}:"$2"_"$3"_"$4; print $0}}' \
        {params.output_prefix}.tmp.legend > {output.legend}
        mv {params.output_prefix}.tmp.hap {output.haps}

        bgzip -f {output.legend}
        bgzip -f {output.haps}
        touch {output.legend}
        touch {output.haps}

        cp {input.sample} {params.output_prefix}.samples

        bcftools convert -H {params.output_prefix} | bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

def get_input_vcfs_as_list(wildcards):
    return(mGen_chunk_vcf_dict[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, mGen_chunk_vcf_dict[str(wildcards.chr)])))

rule merge_1KG_GAMCC_chunks:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz",
    resources: mem = '30G'
    threads: 4
    params: 
        input_string = get_input_vcfs_as_string,
    shell: """
        bcftools concat --ligate-force -Oz -o {output.vcf} {params.input_string}
        tabix {output.vcf}
    """