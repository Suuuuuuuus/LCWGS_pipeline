configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

chromosome = [i for i in range(1,23)]


rule all:
    input:
        loadings1 = "results/population/chip/PCs/loadings.all.tsv",
        loadings2 = "results/population/lc/PCs/loadings.all.tsv",
        # PC4 = "results/population/mixed_all/PCs/PCs.eigenvec",
        # PC1 = "results/population/combined/PCs/PCs.eigenvec",
        # PC3 = "results/population/mixed/PCs/PCs.eigenvec",
        # PC2 = "results/population/chip/PCs/PCs.eigenvec",

rule calculate_PCA_chip:
    input:
        vcf = "results/chip/vcf/chip_qced.vcf.gz",
        vcf2 = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
    output:
        variants = temp("results/population/chip/vcf/variants.tsv"),
        bed = temp("results/population/chip/vcf/chip_pca.bed"),
        bim = temp("results/population/chip/vcf/chip_pca.bim"),
        fam = temp("results/population/chip/vcf/chip_pca.fam"),
        PC = "results/population/chip/PCs/PCs.eigenvec",
        tmp1_vcf = temp("results/population/chip/vcf/tmp1.vcf.gz"),
        tmp2_vcf = temp("results/population/chip/vcf/tmp2.vcf.gz"),
    params:
        num_PCs = 10,
        plink_name = "results/population/chip/vcf/chip_pca",
        PC_name = "results/population/chip/PCs/PCs",
        gam_to_exclude = "data/sample_tsvs/PC_exclude/gam_to_exclude.tsv"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/population/chip/vcf/
        mkdir -p results/population/chip/PCs/

        zgrep -v '^#' {input.vcf2} | cut -f1,2 > {output.variants}

        bcftools view -S ^{params.gam_to_exclude} -Oz -o {output.tmp1_vcf} {input.vcf}
        tabix -f {output.tmp1_vcf} 
        bcftools view -R {output.variants} -Oz -o {output.tmp2_vcf} {output.tmp1_vcf} 
        tabix -f {output.tmp2_vcf} 

        plink --vcf {output.tmp2_vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}

        rm {output.tmp1_vcf}.tbi
        rm {output.tmp2_vcf}.tbi
    """

rule convert_chip_to_bgen:
    input:
        vcf = "results/population/chip/vcf/tmp2.vcf.gz"
    output:
        vcf = temp("results/population/chip/bgen/lc.unphased.vcf.gz"),
        bgen = "results/population/chip/bgen/lc.bgen",
        sample = "results/population/chip/bgen/lc.sample"
    params:
        qctool = tools['qctool']
    localrule: True
    shell: """
        mkdir -p results/population/chip/bgen/

        zcat {input.vcf} | sed 's/|/\//g' | bgzip > {output.vcf}
        
        echo "ID" > {output.sample}
        echo "0" >> {output.sample}
        bcftools query -l {output.vcf} >> {output.sample}

        {params.qctool} -g {output.vcf} -s {output.sample} \
        -og {output.bgen} -bgen-bits 8 -bgen-compression zstd
    """	

rule calculate_loadings_chip:
    input:
        bgen = rules.convert_chip_to_bgen.output.bgen,
        samples = rules.convert_chip_to_bgen.output.sample
    output:
        kinship1 = "results/population/chip/PCs/kinship.all.tsv.gz",
        UDUT1 = "results/population/chip/PCs/UDUT.all.tsv.gz",
        loadings = "results/population/chip/PCs/loadings.all.tsv",
        PCs = "results/population/chip/PCs/PCs.csv",
    params:
        n_PCs = 10,
        qctool = tools['qctool']
    resources:
        mem = '20G'
    threads: 4
    shell: """
        mkdir -p results/population/chip/PCs/

        {params.qctool} \
        -g {input.bgen} -s {input.samples} \
        -kinship {output.kinship1} \
        -UDUT {output.UDUT1} \
        -PCs {params.n_PCs} \
        -osample {output.PCs}

        {params.qctool} \
        -g {input.bgen} -s {input.samples} \
        -load-UDUT {output.UDUT1} \
        -loadings {output.loadings}
    """

rule convert_lc_to_bgen:
    input:
        vcf = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
    output:
        tmp_vcf1 = temp("results/population/lc/bgen/vcf1.gz"),
        tmp_vcf = temp("results/population/lc/bgen/vcf.gz"),
        vcf = temp("results/population/lc/bgen/lc.unphased.vcf.gz"),
        bgen = "results/population/lc/bgen/lc.bgen",
        sample = "results/population/lc/bgen/lc.sample"
    params:
        qctool = tools['qctool'],
        chip_gm_names = "data/sample_tsvs/chip_gm_names.tsv",
        gm_to_exclude = "data/sample_tsvs/PC_exclude/gm_to_exclude.tsv",
    localrule: True
    shell: """
        mkdir -p results/population/lc/bgen/

        bcftools view -S {params.chip_gm_names} -Oz -o {output.tmp_vcf} {input.vcf}
        bcftools view -S ^{params.gm_to_exclude} -Oz -o {output.tmp_vcf1} {output.tmp_vcf}
        tabix -f {output.tmp_vcf1}

        zcat {output.tmp_vcf1} | sed 's/|/\//g' | bgzip > {output.vcf}
        
        echo "ID" > {output.sample}
        echo "0" >> {output.sample}
        bcftools query -l {output.vcf} >> {output.sample}

        {params.qctool} -g {output.vcf} -s {output.sample} \
        -og {output.bgen} -bgen-bits 8 -bgen-compression zstd
    """	

rule calculate_loadings_lc:
    input:
        bgen = rules.convert_lc_to_bgen.output.bgen,
        samples = rules.convert_lc_to_bgen.output.sample
    output:
        kinship1 = "results/population/lc/PCs/kinship.all.tsv.gz",
        UDUT1 = "results/population/lc/PCs/UDUT.all.tsv.gz",
        loadings = "results/population/lc/PCs/loadings.all.tsv",
        PCs = "results/population/lc/PCs/PCs.csv",
    params:
        n_PCs = 10,
        qctool = tools['qctool'],
    resources:
        mem = '20G'
    threads: 4
    shell: """
        mkdir -p results/population/lc/PCs/

        {params.qctool} \
        -g {input.bgen} -s {input.samples} \
        -kinship {output.kinship1} \
        -UDUT {output.UDUT1} \
        -PCs {params.n_PCs} \
        -osample {output.PCs}

        {params.qctool} \
        -g {input.bgen} -s {input.samples} \
        -load-UDUT {output.UDUT1} \
        -loadings {output.loadings}
    """

# rule calculate_PCA_mixed_up:
#     input:
#         vcf1 = "results/chip/vcf/chip_qced.vcf.gz",
#         vcf2 = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
#     output:
#         variants = temp("results/population/mixed/vcf/variants.tsv"),
#         tmp_vcf1 = temp("results/population/mixed/vcf/vcf1.gz"),
#         tmp_vcf2 = temp("results/population/mixed/vcf/vcf2.gz"),
#         tmp_vcf = temp("results/population/mixed/vcf/vcf.gz"),
#         sort_vcf = temp("results/population/mixed/vcf/sorted_vcf.gz"),

#         bed = temp("results/population/mixed/vcf/pca.bed"),
#         bim = temp("results/population/mixed/vcf/pca.bim"),
#         fam = temp("results/population/mixed/vcf/pca.fam"),
#         PC = "results/population/mixed/PCs/PCs.eigenvec",
#     params:
#         chip_samples = "data/sample_tsvs/PC_exclude/chip_half.tsv",
#         lc_samples = "data/sample_tsvs/PC_exclude/lc_half.tsv",
#         tmp_dir = "results/population/mixed/PCs/tmp/",

#         num_PCs = 10,
#         plink_name = "results/population/mixed/vcf/pca",
#         PC_name = "results/population/mixed/PCs/PCs",
#     resources:
#         mem = '10G'
#     shell: """
#         mkdir -p results/population/mixed/PCs/

#         bcftools view -S {params.chip_samples} -Oz -o {output.tmp_vcf1} {input.vcf1}
#         tabix -f {output.tmp_vcf1}
#         bcftools view -S {params.lc_samples} -Oz -o {output.tmp_vcf2} {input.vcf2}
#         tabix -f {output.tmp_vcf2}

#         zgrep -v '^#' {input.vcf2} | cut -f1,2 > {output.variants}
#         bcftools view -R {output.variants} -Oz -o {output.sort_vcf} {output.tmp_vcf1}
#         tabix -f {output.sort_vcf}

#         bcftools merge {output.sort_vcf} {output.tmp_vcf2} |
#         bcftools sort -Oz -o {output.tmp_vcf}
#         tabix -f {output.tmp_vcf}

#         plink --vcf {output.tmp_vcf} --make-bed --out {params.plink_name}
#         plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
#     """

rule calculate_PCA_mixed_up_all:
    input:
        vcf1 = "results/chip/vcf/chip_qced.vcf.gz",
        vcf2 = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
    output:
        variants = temp("results/population/mixed_all/vcf/variants.tsv"),
        tmp_vcf1 = temp("results/population/mixed_all/vcf/vcf1.gz"),
        tmp_vcf2 = temp("results/population/mixed_all/vcf/vcf2.gz"),
        tmp_vcf = temp("results/population/mixed_all/vcf/vcf.gz"),

        bed = temp("results/population/mixed_all/vcf/pca.bed"),
        bim = temp("results/population/mixed_all/vcf/pca.bim"),
        fam = temp("results/population/mixed_all/vcf/pca.fam"),
        PC = "results/population/mixed_all/PCs/PCs.eigenvec",
    params:
        num_PCs = 10,
        plink_name = "results/population/mixed_all/vcf/pca",
        PC_name = "results/population/mixed_all/PCs/PCs",

        chip_gm_names = "data/sample_tsvs/chip_gm_names.tsv",
        gm_to_exclude = "data/sample_tsvs/PC_exclude/gm_to_exclude.tsv",
        gam_to_exclude = "data/sample_tsvs/PC_exclude/gam_to_exclude.tsv"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/population/mixed_all/PCs/
        mkdir -p results/population/mixed_all/vcf/

        zgrep -v '^#' {input.vcf2} | cut -f1,2 > {output.variants}

        bcftools view -S {params.chip_gm_names} -Oz -o {output.tmp_vcf} {input.vcf2}
        bcftools view -S ^{params.gm_to_exclude} -Oz -o {output.tmp_vcf1} {output.tmp_vcf}
        tabix -f {output.tmp_vcf1}

        bcftools view -S ^{params.gam_to_exclude} -Oz -o {output.tmp_vcf} {input.vcf1}
        tabix -f {output.tmp_vcf}
        bcftools view -R {output.variants} -Oz -o {output.tmp_vcf2} {output.tmp_vcf} 
        tabix -f {output.tmp_vcf2} 


        bcftools merge {output.tmp_vcf1} {output.tmp_vcf2} |
        bcftools sort -Oz -o {output.tmp_vcf}
        tabix -f {output.tmp_vcf}

        plink --vcf {output.tmp_vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """