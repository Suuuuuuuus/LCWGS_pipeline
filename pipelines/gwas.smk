configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import io
import os
import re
import json
import pandas as pd
import numpy as np
import math
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]

models = ['lr1']
sources = ['lc', 'chip']

rule concat_chip_sites_vcfs_lc:
    input:
        lc_vcf = expand("results/wip_vcfs/oneKG/vanilla/chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome, allow_missing = True)
    output:
        concat = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
    resources:
        mem = '40G'
    threads: 4
    shell: """
        bcftools concat --threads 4 {input.lc_vcf} | bcftools sort -Oz -o {output.concat}

        gunzip -c {output.concat} | grep '#' > {output.concat}.temp1.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t[%GT\t]\n' \
        {output.concat} >> {output.concat}.temp1.vcf
        bcftools sort -Oz -o {output.concat} {output.concat}.temp1.vcf
        tabix {output.concat}
        rm {output.concat}.temp*
    """

rule calculate_PCA_lc:
    input:
        vcf = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc.vcf.gz"
    output:
        bed = temp("results/wip_vcfs/oneKG/vanilla/chip_sites/lc_pca.bed"),
        bim = temp("results/wip_vcfs/oneKG/vanilla/chip_sites/lc_pca.bim"),
        fam = temp("results/wip_vcfs/oneKG/vanilla/chip_sites/lc_pca.fam"),
        PC = "results/wip_vcfs/oneKG/vanilla/chip_sites/PCs.eigenvec",
        tmp_vcf = temp("results/wip_vcfs/oneKG/vanilla/chip_sites/tmp.vcf.gz")
    params:
        num_PCs = 10,
        plink_name = "results/wip_vcfs/oneKG/vanilla/chip_sites/lc_pca",
        PC_name = "results/wip_vcfs/oneKG/vanilla/chip_sites/PCs",
        rename = "data/rename_tsvs/fv_gm_to_gam.ssv",
        chip_gm_names = "data/sample_tsvs/chip_gm_names.tsv"
    resources:
        mem = '10G'
    shell: """
        bcftools view -S {params.chip_gm_names} {input.vcf} | \
        bcftools reheader -s {params.rename} | \
        bgzip > {output.tmp_vcf}
        tabix {output.tmp_vcf}
        plink --vcf {output.tmp_vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """

rule subset_blood_group_regions_lc:
    input:
        vcf = "results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.dose.vcf.gz"
    output:
        subset_vcf = temp("results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.tmp.vcf.gz"),
        tmp_vcf = temp("results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.tmp1.vcf.gz")
    params:
        regions = "data/blood_group_variants/blood_group_regions.tsv",
        fv_gm_names = "data/sample_tsvs/fv_gm_names.tsv",
        rename = "data/rename_tsvs/fv_gm_to_gam.ssv",
    resources:
        mem = '40G'
    shell: """
        tabix -f {input.vcf}

        bcftools view -R {params.regions} -Oz -o {output.tmp_vcf} {input.vcf}
        tabix -f {output.tmp_vcf}

        bcftools view -S {params.fv_gm_names} {output.tmp_vcf} | \
        bcftools reheader -s {params.rename} | \
        bgzip > {output.subset_vcf}
        tabix -f {output.subset_vcf}
    """

rule prepare_gamcc_samples_lc:
    input:
        PC = "results/wip_vcfs/oneKG/vanilla/chip_sites/PCs.eigenvec"
    output:
        gamcc_gen_samples = "results/gwas/mGenv1_topmed/{model}/lc/vcf/blood_groups.sample"
    params:
        metadata = "data/metadata/GAMCC-lcWGS_metadata.csv",
        chip_gam_names_file = "data/sample_tsvs/chip_gam_names.tsv",
        n = 4
    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[['anon_id', 'Status_v1']]
        chip_gam_names = lcwgsus.read_tsv_as_lst(params.chip_gam_names_file)
        metadata = metadata[metadata['anon_id'].isin(chip_gam_names)]
        metadata.columns = ['ID', 'malaria']
        if wildcards.model == 'mlr':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 1, 'Non-malaria control': 0, 'Severe malaria': 2})
        elif wildcards.model == 'lr1':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 1, 'Non-malaria control': 0, 'Severe malaria': 1})
        elif wildcards.model == 'lr2':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 0, 'Non-malaria control': 0, 'Severe malaria': 1})
        elif wildcards.model == 'lr3':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 2, 'Non-malaria control': 0, 'Severe malaria': 1})
            metadata = metadata[metadata['malaria'] != 2]
        else:
            pass
        metadata = metadata.reset_index(drop = True)

        n_PCs = params.n
        PCs = pd.read_csv(input.PC, sep = ' ', header = None)
        PCs = PCs.iloc[:,1:(n_PCs + 2)]
        PCs.columns = ['ID'] + [f'PC{i}' for i in range(1,(n_PCs + 1))]
        for c in PCs.columns[1:]:
            PCs[c] = PCs[c].astype(str)

        sample = pd.merge(metadata, PCs, on = 'ID')
        if wildcards.model == 'mlr':
            sample.loc[-1] = ['0', 'D', 'C', 'C', 'C', 'C']
        else:
            sample.loc[-1] = ['0', 'B', 'C', 'C', 'C', 'C']
        sample.index = sample.index + 1
        sample = sample.sort_index()
        sample.to_csv(output.gamcc_gen_samples, sep = ' ', index = False, header = True)

rule concat_subset_blood_group_regions_lc:
    input:
        subset_vcf = expand("results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.tmp.vcf.gz", chr = chromosome),
        samples = "results/gwas/mGenv1_topmed/{model}/lc/vcf/blood_groups.sample"
    output:
        concat_vcf = "results/gwas/mGenv1_topmed/{model}/lc/vcf/blood_groups.vcf.gz",
        gen = "results/gwas/mGenv1_topmed/{model}/lc/vcf/blood_groups.gen",
        sample_file = temp("results/gwas/mGenv1_topmed/{model}/lc/vcf/samplefile.tsv"),
        tmp_vcf = temp("results/gwas/mGenv1_topmed/{model}/lc/vcf/blood_groups.tmp.vcf.gz")
    params:
        qctool = tools['qctool']
    shell: """
        mkdir -p results/gwas/mGenv1_topmed/{wildcards.model}/lc/vcf/

        bcftools concat {input.subset_vcf} | bcftools sort -Oz -o {output.concat_vcf}
        tabix -f {output.concat_vcf}

        cat {input.samples} | tail -n +3 | cut -d ' ' -f1 > {output.sample_file}
        bcftools view -S {output.sample_file} -Oz -o {output.tmp_vcf} {output.concat_vcf}
        tabix -f {output.tmp_vcf}

        {params.qctool} \
        -g {output.tmp_vcf} -s {input.samples} -vcf-genotype-field GP -og {output.gen}
    """

rule blood_group_gwas:
    input:
        gen = "results/gwas/mGenv1_topmed/{model}/{source}/vcf/blood_groups.gen",
        samples = "results/gwas/mGenv1_topmed/{model}/{source}/vcf/blood_groups.sample"
    output:
        result = "results/gwas/mGenv1_topmed/{model}/{source}/results/stats.out"
    params:
        snptest = tools['snptest'],
        pheno = "malaria"
    localrule: True
    shell: """
        mkdir -p results/gwas/mGenv1_topmed/{wildcards.model}/{wildcards.source}/results/

        if [[ {wildcards.model} == "mlr" ]]
        then
            {params.snptest} \
            -data {input.gen} {input.samples} \
            -o {output.result} \
            -frequentist dom \
            -method newml \
            -pheno {params.pheno} \
            -baseline_phenotype 0
        else
            {params.snptest} \
            -data {input.gen} {input.samples} \
            -o {output.result} \
            -frequentist dom \
            -method newml \
            -pheno {params.pheno}
        fi
    """


### Below for chip

rule calculate_PCA_chip:
    input:
        vcf = "results/chip/vcf/chip_qced.vcf.gz"
    output:
        bed = temp("results/chip/PCs/chip_pca.bed"),
        bim = temp("results/chip/PCs/chip_pca.bim"),
        fam = temp("results/chip/PCs/chip_pca.fam"),
        PC = "results/chip/PCs/PCs.eigenvec"
    params:
        num_PCs = 10,
        plink_name = "results/chip/PCs/chip_pca",
        PC_name = "results/chip/PCs/PCs"
    resources:
        mem = '10G'
    localrule: True
    shell: """
        plink --vcf {input.vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """

rule subset_blood_group_regions_chip:
    input:
        vcf = "results/chip/imputed/topmed/vcf/chr{chr}.dose.vcf.gz"
    output:
        subset_vcf = temp("results/chip/imputed/topmed/vcf/chr{chr}.tmp.vcf.gz")
    params:
        regions = "data/blood_group_variants/blood_group_regions.tsv"
    resources:
        mem = '40G'
    shell: """
        tabix -f {input.vcf}

        bcftools view -R {params.regions} -Oz -o {output.subset_vcf} {input.vcf}
        tabix -f {output.subset_vcf}
    """

rule prepare_gamcc_samples_chip:
    input:
        PC = "results/chip/PCs/PCs.eigenvec"
    output:
        gamcc_gen_samples = "results/gwas/mGenv1_topmed/{model}/chip/vcf/blood_groups.sample"
    params:
        metadata = "data/metadata/GAMCC-lcWGS_metadata.csv",
        chip_gam_names_file = "data/sample_tsvs/chip_gam_names.tsv",
        n = 4
    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[['anon_id', 'Status_v1']]
        chip_gam_names = lcwgsus.read_tsv_as_lst(params.chip_gam_names_file)
        metadata = metadata[metadata['anon_id'].isin(chip_gam_names)]
        metadata.columns = ['ID', 'malaria']
        if wildcards.model == 'mlr':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 1, 'Non-malaria control': 0, 'Severe malaria': 2})
        elif wildcards.model == 'lr1':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 1, 'Non-malaria control': 0, 'Severe malaria': 1})
        elif wildcards.model == 'lr2':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 0, 'Non-malaria control': 0, 'Severe malaria': 1})
        elif wildcards.model == 'lr3':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 2, 'Non-malaria control': 0, 'Severe malaria': 1})
            metadata = metadata[metadata['malaria'] != 2]
        else:
            pass
        metadata = metadata.reset_index(drop = True)

        n_PCs = params.n
        PCs = pd.read_csv(input.PC, sep = ' ', header = None)
        PCs = PCs.iloc[:,1:(n_PCs + 2)]
        PCs.columns = ['ID'] + [f'PC{i}' for i in range(1,(n_PCs + 1))]
        for c in PCs.columns[1:]:
            PCs[c] = PCs[c].astype(str)

        sample = pd.merge(metadata, PCs, on = 'ID')
        if wildcards.model == 'mlr':
            sample.loc[-1] = ['0', 'D', 'C', 'C', 'C', 'C']
        else:
            sample.loc[-1] = ['0', 'B', 'C', 'C', 'C', 'C']
        sample.index = sample.index + 1
        sample = sample.sort_index()
        sample.to_csv(output.gamcc_gen_samples, sep = ' ', index = False, header = True)

rule concat_subset_blood_group_regions_chip:
    input:
        subset_vcf = expand("results/chip/imputed/topmed/vcf/chr{chr}.tmp.vcf.gz", chr = chromosome),
        samples = "results/gwas/mGenv1_topmed/{model}/chip/vcf/blood_groups.sample"
    output:
        concat_vcf = "results/gwas/mGenv1_topmed/{model}/chip/vcf/blood_groups.vcf.gz",
        gen = "results/gwas/mGenv1_topmed/{model}/chip/vcf/blood_groups.gen",
        sample_file = temp("results/gwas/mGenv1_topmed/{model}/chip/vcf/samplefile.tsv"),
        tmp_vcf = temp("results/gwas/mGenv1_topmed/{model}/chip/vcf/blood_groups.tmp.vcf.gz")
    params:
        qctool = tools['qctool']
    shell: """
        mkdir -p results/gwas/mGenv1_topmed/{wildcards.model}/chip/vcf/

        bcftools concat {input.subset_vcf} | bcftools sort -Oz -o {output.concat_vcf}
        tabix -f {output.concat_vcf}

        cat {input.samples} | tail -n +3 | cut -d ' ' -f1 > {output.sample_file}
        bcftools view -S {output.sample_file} -Oz -o {output.tmp_vcf} {output.concat_vcf}
        tabix -f {output.tmp_vcf}

        {params.qctool} \
        -g {output.tmp_vcf} -s {input.samples} -vcf-genotype-field GP -og {output.gen}
    """