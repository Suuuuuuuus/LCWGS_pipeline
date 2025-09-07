configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

chromosome = [i for i in range(1,23)]
imputation_dir_new = config['imputation_dir'][:2]

rule all:
    input:
        PC = "results/population/combined/PCs/PCs.eigenvec"

rule concat_onekg_vcf:
    input:
        onekg_vcf = expand("data/ref_panel/oneKG_30x/oneKG.chr{chr}.vcf.gz", chr = chromosome),
    output:
        concat_vcf = temp("data/ref_panel/oneKG_30x/oneKG.vcf.gz")
    resources:
        mem = '40G'
    threads: 4
    params:
        samples = '/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/samples_mandinka.tsv'
    shell: """
        bcftools concat {input.onekg_vcf} | \
        bcftools view -S {params.samples} | \
        bcftools view -i 'MAF>0.05' | \
        bcftools sort -Oz -o {output.concat_vcf}
        tabix -f {output.concat_vcf}
    """

rule concat_imputed_vcf:
    input:
        onekg_vcf = expand("results/wip_vcfs/oneKG/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz", chr = chromosome),
    output:
        concat_vcf = temp("results/wip_vcfs/oneKG/vanilla/high_info_high_af/lc.vcf.gz")
    resources:
        mem = '40G'
    threads: 4
    shell: """
        bcftools concat {input.onekg_vcf} | \
        bcftools view -i 'MAF>0.05' | \
        bcftools sort -Oz -o {output.concat_vcf}
        tabix -f {output.concat_vcf}
    """

rule merge_batch_vcf_with_1kg:
    input:
        onekg_vcf = "data/ref_panel/oneKG_30x/oneKG.vcf.gz",
        merged = "results/wip_vcfs/oneKG/vanilla/high_info_high_af/lc.vcf.gz"
    output:
        vcf = "results/population/combined/vcf/all.vcf.gz"
    resources:
        mem = '40G'
    threads: 4
    params: 
        tmp_dir = "results/population/combined/PCs/tmp/"
    shell: """
        mkdir -p {params.tmp_dir}
        mkdir -p results/population/combined/vcf/

        bcftools isec -n =2 {input.onekg_vcf} {input.merged} -w1 -Oz -p {params.tmp_dir}
        bcftools merge {input.onekg_vcf} {input.merged} --regions-file {params.tmp_dir}sites.txt | \
        bcftools annotate -x ID -I +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o {output.vcf}

        tabix -f {output.vcf}
        rm -r {params.tmp_dir}
    """

rule calculate_PCA:
    input:
        vcf = "results/population/combined/vcf/all.vcf.gz"
    output:
        bed = temp("results/population/combined/vcf/pca.bed"),
        bim = temp("results/population/combined/vcf/pca.bim"),
        fam = temp("results/population/combined/vcf/pca.fam"),
        PC = "results/population/combined/PCs/PCs.eigenvec"
    params:
        num_PCs = 10,
        plink_name = "results/population/combined/vcf/pca",
        PC_name = "results/population/combined/PCs/PCs"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/population/combined/PCs/

        plink --vcf {input.vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """