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

# The followings are global parameters from `activate`:
QUILT_HOME = config["QUILT_HOME"]
ANALYSIS_DIR = config["ANALYSIS_DIR"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

rule calculate_imputation_accuracy:
    input:
        quilt_vcf = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz",
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "results/imputation_metrics/lc_chip/by_variant/lc.chip.typed.chr{chr}.h.tsv",
        h_impacc = "results/imputation_metrics/lc_chip/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv",
        v_report = "results/imputation_metrics/lc_chip/by_sample/lc.chip.typed.chr{chr}.v.tsv",
        v_impacc = "results/imputation_metrics/lc_chip/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker']
    run:
        mini = False
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'
        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af
        sample_linker = pd.read_table(params.linker, sep = ',')
        # by_variant = True # True if average each row, False if average each column (sample)
        
        af = lcwgsus.read_af(af)

        MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.95, 1])
        
        chromosome = wildcards.chr
        common_cols = ['chr', 'pos', 'ref', 'alt']

        if not mini:
            sample_linker = sample_linker[~sample_linker['Sample_Name'].str.contains('mini')]
        else:
            sample_linker = sample_linker[sample_linker['Sample_Name'].str.contains('mini')]
        rename_map = dict(zip(sample_linker['Sample_Name'], sample_linker['Chip_Name']))

        lc = lcwgsus.read_vcf(quilt_vcf).sort_values(by = ['chr', 'pos'])
        chip = lcwgsus.read_vcf(chip_vcf).sort_values(by = ['chr', 'pos'])

        if not mini:
            lc = lc.drop(columns = lc.columns[lc.columns.str.contains('mini')])
        else:
            lc = lc.drop(columns = lc.columns[~lc.columns.str.contains('mini')])

        lc = lc.drop(columns = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
        chip = chip.drop(columns = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

        res = lcwgsus.intersect_dfs([chip, lc, af])
        chip = res[0]
        lc = res[1]
        af = res[2]

        lc_samples = lc.columns[lc.columns.str.contains(lc_sample_prefix)]
        chip_samples = chip.columns[chip.columns.str.contains(chip_sample_prefix)]
        lc_to_retain = lcwgsus.find_matching_samples(lc_samples, chip_samples, rename_map)
        quilt_columns = list(lc.columns[~lc.columns.str.contains(lc_sample_prefix)]) + lc_to_retain
        lc = lc[quilt_columns]

        lc = lc.apply(lcwgsus.extract_DS, axis = 1)
        chip = chip.apply(lcwgsus.encode_genotype, axis = 1)

        chip_order = []
        for i in lc_to_retain:
            chip_order.append(rename_map[i])
        chip = chip[common_cols + chip_order]
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = 'results/imputation_metrics/lc_chip/by_variant/', 
                                                   save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/by_variant/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

'''
rule plot_imputation_accuracy:
    input:
        r2 = expand("results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv", id = seq_to_extract, panel = panels)
    output:
        graph = "graphs/imputation_vs_chip.png"
    resources:
        mem = '10G'
    params:
        samples = seq_to_extract,
        panels = panels,
        linker = config['sample_linker']
    run:
        samples = params.samples
        panels = params.panels
        linker = pd.read_table(params.linker, sep = ',')
        fv = linker[(linker['Seq_Name'].isin(samples)) & (~linker['Sample_Name'].str.contains('mini'))]
        mini = linker[(linker['Seq_Name'].isin(samples)) & (linker['Sample_Name'].str.contains('mini'))]
        samples_fv = fv['Seq_Name'].to_list()
        samples_mini = mini['Seq_Name'].to_list()
        r2_fv, _ = lcwgsus.read_r2(panels, samples_fv)
        r2_fv = lcwgsus.aggregate_r2(r2_fv)
        r2_mini, _ = lcwgsus.read_r2(panels, samples_mini)
        r2_mini = lcwgsus.aggregate_r2(r2_mini)

        plt.figure(figsize = (10,6))
        for df in r2_fv:
            panel = df['panel'].values[0]
            plt.plot(np.arange(1, df.shape[0]+1), df['corr'], label = panel)
        for df in r2_mini:
            panel = df['panel'].values[0]
            plt.plot(np.arange(1, df.shape[0]+1), df['corr'], label = panel + '_mini', ls ='--')
        plt.xticks(np.arange(1, r2_fv[0].shape[0]+1), r2_fv[0]['AF'], rotation = 45)
        plt.xlabel('Allele frequencies (%)')
        plt.legend()
        plt.text(x = -1.5, y = 1.04, s = 'Aggregated imputation accuracy ($r^2$)')
        plt.grid(alpha = 0.5)
        plt.savefig(output.graph, bbox_inches = "tight", dpi=300)
'''
