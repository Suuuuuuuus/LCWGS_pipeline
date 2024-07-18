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
from lcwgsus.variables import *

sample_linker = pd.read_table(config['sample_linker'], sep = ',')
chromosome = [i for i in range(1,23)]

imp_dir = config['imputation_dir']
case_controls = ['non-malaria_control', 'mild_malaria', 'severe_malaria']
ethnicities = ['fula', 'jola', 'mandinka', 'wollof']
pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['imputation_dir']
lc_vcf_dir = config['lc_vcf_dir']
hc_vcf_dir = config['hc_vcf_dir']

def get_lc_vcf_dir(wildcards):
    ix = imputation_dir.index(wildcards.imp_dir)
    return lc_vcf_dir[ix]

def get_hc_vcf_dir(wildcards):
    ix = imputation_dir.index(wildcards.imp_dir)
    return hc_vcf_dir[ix]

rule copy_vcf_in_working_dir:
    input:
        lc_vcf_dir = get_lc_vcf_dir,
        hc_vcf_dir = get_hc_vcf_dir
    output:
        vcfs = temp(expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, allow_missing = True))
    localrule: True
    shell: """
        mkdir -p {wildcards.imp_dir}vcf/
        mkdir -p {wildcards.imp_dir}impacc/
        mkdir -p {wildcards.imp_dir}vcf/all_samples/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/all_samples/hc_vcf/

        for i in {{1..22}}
        do
            cp -f {input.lc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/lc_vcf/lc.chr$i.vcf.gz
            cp -f {input.hc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/hc_vcf/hc.chr$i.vcf.gz
        done
    """

rule split_vcf_by_eth:
    input:
        vcf = '{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz'
    output:
        eth_vcf = '{imp_dir}vcf/by_eth/{pair}_vcf/{eth}.{pair}.chr{chr}.vcf.gz'
    resources:
        mem = '10G'
    shell: """
        mkdir -p {wildcards.imp_dir}vcf/by_eth/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_eth/hc_vcf/

        bcftools view -S data/file_lsts/samples_subset/by_ethnicity/{wildcards.eth}_samples_{wildcards.pair}.tsv \
        -Oz -o {output.eth_vcf} {input.vcf}
    """

rule split_vcf_by_cc:
    input:
        vcf = '{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz'
    output:
        cc_vcf = '{imp_dir}vcf/by_cc/{pair}_vcf/{cc}.{pair}.chr{chr}.vcf.gz'
    resources:
        mem = '10G'
    shell: """
        mkdir -p {wildcards.imp_dir}vcf/by_cc/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_cc/hc_vcf/
        
        bcftools view -S data/file_lsts/samples_subset/by_case_control/{wildcards.cc}_samples_{wildcards.pair}.tsv \
        -Oz -o {output.cc_vcf} {input.vcf}
    """

rule subset_lc_samples:
    input:
        quilt_vcf = '{imp_dir}vcf/all_samples/lc_vcf/lc.chr{chr}.vcf.gz',
        chip_vcf = '{imp_dir}vcf/all_samples/hc_vcf/hc.chr{chr}.vcf.gz'
    output:
        ss_lc_vcf = temp('{imp_dir}vcf/all_samples/lc_vcf/lc.subset.chr{chr}.vcf.gz'),
        ss_hc_vcf = temp('{imp_dir}vcf/all_samples/hc_vcf/hc.subset.chr{chr}.vcf.gz'),
        tmp_names = temp('{imp_dir}vcf/all_samples/samples.chr{chr}.tsv')
    resources:
        mem = '10G'
    run: 
        hc_names = lcwgsus.bcftools_get_samples(input.chip_vcf)
        lc_names = lcwgsus.bcftools_get_samples(input.quilt_vcf)

        rename_map = lcwgsus.generate_rename_map()
        lc = wildcards.imp_dir.split('/')[-2].split('_')[0]
        samples = lcwgsus.find_matching_samples(hc_names, rename_map, lc = lc)
        lcwgsus.save_lst(output.tmp_names, samples)

        shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_lc_vcf} {input.quilt_vcf}")
        shell("cp {input.chip_vcf} {output.ss_hc_vcf}")

rule calculate_imputation_accuracy_all:
    input:
        quilt_vcf = rules.subset_lc_samples.output.ss_lc_vcf,
        chip_vcf = rules.subset_lc_samples.output.ss_hc_vcf,
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv",
        h_impacc = "{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv",
        v_report = "{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv",
        v_impacc = "{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv",
        af = "{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv"
    resources:
        mem = '100G'
    threads: 16
    params:
        linker = config['sample_linker'],
        common_outdir = "{imp_dir}impacc/all_samples/",
        common_savedir = "{imp_dir}vcf/all_samples/",
        save_filtered_files = False,
        lc_vcf = "{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz",
        hc_vcf = "{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz",
        chrom = "{chr}"
    run:
        mini = False

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af_txt = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt, save_vcfs = params.save_filtered_files, lc_vcf_outdir = params.common_savedir + "filtered_vcfs/", hc_vcf_outdir = params.common_savedir + "filtered_vcfs/", lc_vcf_name = "lc.chr" + wildcards.chr + ".vcf.gz", hc_vcf_name = "hc.chr" + wildcards.chr + ".vcf.gz", af_outdir = params.common_savedir + "af/", af_name = "af.chr" + wildcards.chr + ".tsv")

        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = 'chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = COMMON_COLS)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = 'chr' + wildcards.chr +'.h.impacc.tsv')

        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report,
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable
        if params.save_filtered_files:
            lcwgsus.rezip_vcf(params.lc_vcf)
            lcwgsus.rezip_vcf(params.hc_vcf)

rule get_badly_imputed_variants:
    input:
        h_impaccs = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        tsv = "{imp_dir}impacc/all_samples/by_variant/all_r2less0.5.tsv"
    resources:
        mem = '30G'
    params:
        threshold = 0.5
    run:
        h_lst = input.h_impaccs
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        lcwgsus.get_badly_imputed_regions(
            wildcards.imp_dir + 'impacc/all_samples/by_variant/',
            retain_cols = COMMON_COLS + ['MAF', 'r2'],
            threshold = params.threshold,
            save_file = True,
            outdir = wildcards.imp_dir + 'impacc/all_samples/by_variant/',
            save_name = 'all_r2less0.5.tsv'
        )

rule aggregate_impaccs_all:
    input:
        h_impaccs = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        h_impacc = "{imp_dir}impacc/all_samples/by_variant/all.h.impacc.tsv",
        v_impacc = "{imp_dir}impacc/all_samples/by_sample/all.v.impacc.tsv"
    resources:
        mem = '30G'
    params:
        impacc_outdir = "{imp_dir}impacc/all_samples/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.impacc_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs, save_file = True,
        outdir = outdir_h, save_name = 'all.h.impacc.tsv')

        v_lst = input.v_impaccs
        outdir_v = params.impacc_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs, save_file = True,
        outdir = outdir_v, save_name = 'all.v.impacc.tsv')

rule plot_imputation_accuracy_all:
    input:
        h_impacc = rules.aggregate_impaccs_all.output.h_impacc,
        v_impacc = rules.aggregate_impaccs_all.output.v_impacc
    output:
        r2NRC_h = "{imp_dir}graphs/all_samples/by_variant/r2_NRC.png",
        ccd_h = "{imp_dir}graphs/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = "{imp_dir}graphs/all_samples/by_sample/r2_NRC.png",
        ccd_v = "{imp_dir}graphs/all_samples/by_sample/ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = "{imp_dir}graphs/all_samples/"
    run:
        h = pd.read_csv(input.h_impacc, sep = '\t')
        outdir_h = params.common_outdir + "by_variant/"

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = "r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = "ccd_by_genotype.png")

        v = pd.read_csv(input.v_impacc, sep = '\t')
        outdir_v = params.common_outdir + "by_sample/"

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = "r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = "ccd_by_genotype.png")

rule calculate_imputation_sumstat:
    input:
        h_impaccs = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        sumstats = "{imp_dir}summary_metrics.tsv"
    resources:
        mem = '30G'
    params:
        axis = "v"
    run:
        lcwgsus.calculate_imputation_sumstats(wildcards.imp_dir, subset = False, axis = params.axis, save_file = True)

rule subset_lc_samples_by_eth:
    input:
        quilt_vcf = '{imp_dir}vcf/by_eth/lc_vcf/{eth}.lc.chr{chr}.vcf.gz',
        chip_vcf = '{imp_dir}vcf/by_eth/hc_vcf/{eth}.hc.chr{chr}.vcf.gz'
    output:
        ss_lc_vcf = temp('{imp_dir}vcf/by_eth/hc_vcf/{eth}.lc.subset.chr{chr}.vcf.gz'),
        tmp_names = temp('{imp_dir}vcf/by_eth/{eth}.samples.chr{chr}.tsv')
    resources:
        mem = '10G'
    run: 
        lc_names = lcwgsus.bcftools_get_samples(input.quilt_vcf)
        chip_names = lcwgsus.bcftools_get_samples(input.chip_vcf)
        sample_linker = pd.read_csv(config['sample_linker'])
        sample_linker = sample_linker[sample_linker['Sample_Name'].isin(lc_names)]
        d = {k:v for k, v in zip(sample_linker['Chip_Name'].values, sample_linker['Sample_Name'].values)}
        lc_return = []
        for i in chip_names:
            lc_return.append(d[i])

        lcwgsus.save_lst(output.tmp_names, lc_return)
        shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_lc_vcf} {input.quilt_vcf}")

rule calculate_imputation_accuracy_by_eth:
    input:
        quilt_vcf = rules.subset_lc_samples_by_eth.output.ss_lc_vcf,
        chip_vcf = '{imp_dir}vcf/by_eth/hc_vcf/{eth}.hc.chr{chr}.vcf.gz',
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.tsv",
        h_impacc = "{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv",
        v_report = "{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.tsv",
        v_impacc = "{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker'],
        common_outdir = "{imp_dir}impacc/by_eth/"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af_txt = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = wildcards.eth + '.chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule aggregate_impaccs_by_eth:
    input:
        h_impaccs = expand("{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        h_impacc = "{imp_dir}impacc/by_eth/by_variant/{eth}.all.h.impacc.tsv",
        v_impacc = "{imp_dir}impacc/by_eth/by_sample/{eth}.all.v.impacc.tsv"
    resources:
        mem = '30G'
    params:
        impacc_outdir = "{imp_dir}impacc/by_eth/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.impacc_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs, save_file = True,
        outdir = params.impacc_outdir + 'by_variant/', save_name = wildcards.eth + '.all.h.impacc.tsv')

        v_lst = input.v_impaccs
        outdir_v = params.impacc_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs, save_file = True,
        outdir = params.impacc_outdir + 'by_sample/', save_name = wildcards.eth + '.all.v.impacc.tsv')

rule plot_imputation_accuracy_by_eth:
    input:
        h_impacc = rules.aggregate_impaccs_by_eth.output.h_impacc,
        v_impacc = rules.aggregate_impaccs_by_eth.output.v_impacc
    output:
        r2NRC_h = "{imp_dir}graphs/by_eth/by_variant/{eth}.r2_NRC.png",
        ccd_h = "{imp_dir}graphs/by_eth/by_variant/{eth}.ccd_by_genotype.png",
        r2NRC_v = "{imp_dir}graphs/by_eth/by_sample/{eth}.r2_NRC.png",
        ccd_v = "{imp_dir}graphs/by_eth/by_sample/{eth}.ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = "{imp_dir}graphs/by_eth/"
    run:
        h = pd.read_csv(input.h_impacc, sep = '\t')
        outdir_h = params.common_outdir + "by_variant/"

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = wildcards.eth + ".r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = wildcards.eth + ".ccd_by_genotype.png")

        v = pd.read_csv(input.v_impacc, sep = '\t')
        outdir_v = params.common_outdir + "by_sample/"

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = wildcards.eth + ".r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = wildcards.eth + ".ccd_by_genotype.png")

rule subset_lc_samples_by_cc:
    input:
        quilt_vcf = '{imp_dir}vcf/by_cc/lc_vcf/{cc}.lc.chr{chr}.vcf.gz',
        chip_vcf = '{imp_dir}vcf/by_cc/hc_vcf/{cc}.hc.chr{chr}.vcf.gz'
    output:
        ss_lc_vcf = temp('{imp_dir}vcf/by_cc/hc_vcf/{cc}.lc.subset.chr{chr}.vcf.gz'),
        tmp_names = temp('{imp_dir}vcf/by_cc/{cc}.samples.chr{chr}.tsv')
    resources:
        mem = '10G'
    run: 
        lc_names = lcwgsus.bcftools_get_samples(input.quilt_vcf)
        chip_names = lcwgsus.bcftools_get_samples(input.chip_vcf)
        sample_linker = pd.read_csv(config['sample_linker'])
        sample_linker = sample_linker[sample_linker['Sample_Name'].isin(lc_names)]
        d = {k:v for k, v in zip(sample_linker['Chip_Name'].values, sample_linker['Sample_Name'].values)}
        lc_return = []
        for i in chip_names:
            lc_return.append(d[i])

        lcwgsus.save_lst(output.tmp_names, lc_return)
        shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_lc_vcf} {input.quilt_vcf}")

rule calculate_imputation_accuracy_by_cc:
    input:
        quilt_vcf = rules.subset_lc_samples_by_cc.output.ss_lc_vcf,
        chip_vcf = '{imp_dir}vcf/by_cc/hc_vcf/{cc}.hc.chr{chr}.vcf.gz',
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.tsv",
        h_impacc = "{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv",
        v_report = "{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.tsv",
        v_impacc = "{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker'],
        common_outdir = "{imp_dir}impacc/by_cc/"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af_txt = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = wildcards.cc + '.chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule aggregate_impaccs_by_cc:
    input:
        h_impaccs = expand("{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        h_impacc = "{imp_dir}impacc/by_cc/by_variant/{cc}.all.h.impacc.tsv",
        v_impacc = "{imp_dir}impacc/by_cc/by_sample/{cc}.all.v.impacc.tsv"
    resources:
        mem = '30G'
    params:
        impacc_outdir = "{imp_dir}impacc/by_cc/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.impacc_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs, save_file = True,
        outdir = params.impacc_outdir + 'by_variant/', save_name = wildcards.cc + '.all.h.impacc.tsv')

        v_lst = input.v_impaccs
        outdir_v = params.impacc_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs, save_file = True,
        outdir = params.impacc_outdir + 'by_sample/', save_name = wildcards.cc + '.all.v.impacc.tsv')

rule plot_imputation_accuracy_by_cc:
    input:
        h_impacc = rules.aggregate_impaccs_by_cc.output.h_impacc,
        v_impacc = rules.aggregate_impaccs_by_cc.output.v_impacc
    output:
        r2NRC_h = "{imp_dir}graphs/by_cc/by_variant/{cc}.r2_NRC.png",
        ccd_h = "{imp_dir}graphs/by_cc/by_variant/{cc}.ccd_by_genotype.png",
        r2NRC_v = "{imp_dir}graphs/by_cc/by_sample/{cc}.r2_NRC.png",
        ccd_v = "{imp_dir}graphs/by_cc/by_sample/{cc}.ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = "{imp_dir}graphs/by_cc/"
    run:
        h = pd.read_csv(input.h_impacc, sep = '\t')
        outdir_h = params.common_outdir + "by_variant/"

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = wildcards.cc + ".r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = wildcards.cc + ".ccd_by_genotype.png")

        v = pd.read_csv(input.v_impacc, sep = '\t')
        outdir_v = params.common_outdir + "by_sample/"

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = wildcards.cc + ".r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = wildcards.cc + ".ccd_by_genotype.png")

rule calculate_imputation_sumstat_all:
    input:
        cc_h_impaccs = expand("{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv", chr = chromosome, cc = case_controls, allow_missing = True),
        cc_v_impaccs = expand("{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv", chr = chromosome, cc = case_controls, allow_missing = True),
        eth_h_impaccs = expand("{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv", chr = chromosome, eth = ethnicities, allow_missing = True),
        eth_impaccs = expand("{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv", chr = chromosome, eth = ethnicities, allow_missing = True)
    output:
        sumstats = "{imp_dir}summary_metrics_all.tsv"
    resources:
        mem = '30G'
    params:
        axis = "v"
    run:
        lcwgsus.calculate_imputation_sumstats(wildcards.imp_dir, subset = True, axis = params.axis, save_file = True, save_name = "summary_metrics_all.tsv")
