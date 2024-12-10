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
from hla_phase import *
from hla_align_functions import *
from hla_align import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
v3570_db_alleles = [5093, 6106, 5657, 869, 624]

rule concat_chip_sites_vcfs:
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

rule calculate_PCA:
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
        fv_gm_names = "data/sample_tsvs/fv_gm_names.tsv"
    resources:
        mem = '10G'
    shell: """
        bcftools view -S {params.fv_gm_names} {input.vcf} | \
        bcftools reheader -s {params.rename} | \
        bgzip > {output.tmp_vcf}
        tabix {output.tmp_vcf}
        plink --vcf {output.tmp_vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """

rule subset_blood_group_regions:
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

rule prepare_gamcc_samples:
    input:
        PC = "results/wip_vcfs/oneKG/vanilla/chip_sites/PCs.eigenvec"
    output:
        gamcc_gen_samples = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.sample"
    params:
        metadata = "data/metadata/GAMCC-lcWGS_metadata.csv",
        fv_gam_names_file = "data/sample_tsvs/fv_gam_names.tsv",
        n = 4
    run:
        metadata = pd.read_csv(params.metadata)
        metadata = metadata[['anon_id', 'Status_v1']]
        fv_gam_names = lcwgsus.read_tsv_as_lst(params.fv_gam_names_file)
        metadata = metadata[metadata['anon_id'].isin(fv_gam_names)]
        metadata.columns = ['ID', 'malaria']
        if wildcards.model == 'mlr':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 'mild', 'Non-malaria control': 'control', 'Severe malaria': 'severe'})
        elif wildcards.model == 'lr1':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 'case', 'Non-malaria control': 'control', 'Severe malaria': 'case'})
        elif wildcards.model == 'lr2':
            metadata['malaria'] = metadata['malaria'].map({'Mild malaria': 'control', 'Non-malaria control': 'control', 'Severe malaria': 'case'})
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

rule concat_subset_blood_group_regions:
    input:
        subset_vcf = expand("results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.tmp.vcf.gz", chr = chromosome),
        samples = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.sample"
    output:
        concat_vcf = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.vcf.gz",
        gen = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.gen",
        sample_file = temp("results/gwas/mGenv1_topmed/{model}/vcf/samplefile.tsv"),
        tmp_vcf = temp("results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.tmp.vcf.gz")
    params:
        qctool = tools['qctool']
    shell: """
        mkdir -p results/gwas/mGenv1_topmed/{wildcards.model}/vcf/

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
        gen = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.gen",
        samples = "results/gwas/mGenv1_topmed/{model}/vcf/blood_groups.sample"
    output:
        result = "results/gwas/mGenv1_topmed/{model}/results/stats.out"
    params:
        snptest = tools['snptest'],
        pheno = "malaria"
    shell: """
        mkdir -p results/gwas/mGenv1_topmed/{wildcards.model}/results/

        {params.snptest} \
        -data {input.gen} {input.samples} \
        -o {output.result} \
        -frequentist gen \
        -method newml \
        -pheno {params.pheno} \
        -baseline_phenotype control
    """