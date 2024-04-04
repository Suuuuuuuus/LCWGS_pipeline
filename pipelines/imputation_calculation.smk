include: "auxiliary.smk"
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

rule calculate_imputation_accuracy:
    input:
        quilt_vcf = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz",
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.tsv",
        h_impacc = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv",
        v_report = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.tsv",
        v_impacc = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker']
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                                   save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy:
    input:
        h_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv", chr = chromosome),
        v_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv", chr = chromosome)
    output:
        r2NRC_h = "graphs/imputation/lc_chip/all_samples/by_variant/r2_NRC.png",
        ccd_h = "graphs/imputation/lc_chip/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = "graphs/imputation/lc_chip/all_samples/by_sample/r2_NRC.png",
        ccd_v = "graphs/imputation/lc_chip/all_samples/by_sample/ccd_by_genotype.png"
    resources:
        mem = '30G'
    run:
        h_lst = input.h_impaccs
        v_lst = input.v_impaccs
        outdir_h = "graphs/imputation/lc_chip/all_samples/by_variant/"
        outdir_v = "graphs/imputation/lc_chip/all_samples/by_sample/"

        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = "r2_NRC.png")

        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = "ccd_by_genotype.png")
    
        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = "r2_NRC.png")

        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = "ccd_by_genotype.png")

rule subset_vcfs:
    input:
        quilt_vcf = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz"
    output:
        eth_quilt_vcfs = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        eth_chip_vcfs = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz",

rule calculate_imputation_accuracy_by_ethnic:
    input:
        quilt_vcf = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz",
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.tsv",
        h_impacc = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv",
        v_report = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.tsv",
        v_impacc = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker']
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                                   save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy:
    input:
        h_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv", chr = chromosome),
        v_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv", chr = chromosome)
    output:
        r2NRC_h = "graphs/imputation/lc_chip/all_samples/by_variant/r2_NRC.png",
        ccd_h = "graphs/imputation/lc_chip/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = "graphs/imputation/lc_chip/all_samples/by_sample/r2_NRC.png",
        ccd_v = "graphs/imputation/lc_chip/all_samples/by_sample/ccd_by_genotype.png"
    resources:
        mem = '30G'
    run:
        h_lst = input.h_impaccs
        v_lst = input.v_impaccs
        outdir_h = "graphs/imputation/lc_chip/all_samples/by_variant/"
        outdir_v = "graphs/imputation/lc_chip/all_samples/by_sample/"

        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = "r2_NRC.png")

        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = "ccd_by_genotype.png")
    
        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = "r2_NRC.png")

        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = "ccd_by_genotype.png")

