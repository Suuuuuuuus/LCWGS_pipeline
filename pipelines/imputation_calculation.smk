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

rule get_chip_vcf:
    input:
        chip_qced = "results/chip/vcf/chip_qced.vcf.gz"
    output:
        chip_vcf = temp("results/chip/tmp/{id}/{id}.vcf.gz")
    params:
        seq_name = lambda wildcards: wildcards.id,
        chip_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Chip_Name'].values[0]
    resources:
        mem_mb = 30000
    shell: """
        mkdir -p results/chip/tmp/{wildcards.id}/

        if (bcftools query -l {input.chip_qced} | grep -q {params.chip_name}); then
            bcftools view -s {params.chip_name} -Oz -o {output.chip_vcf} {input.chip_qced}
        else
            touch {output.chip_vcf}
        fi
    """

rule get_imputation_vcf:
    input:
        imputation_result = "results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz"
    output:
        imputation_vcf = temp("results/imputation/tmp/{id}/{panel}_chr{chr}.vcf.gz")
    params:
        seq_name = lambda wildcards: wildcards.id,
        sample_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Sample_Name'].values[0]
    resources:
        mem_mb = 30000
    shell: """
        if (bcftools query -l {input.imputation_result} | grep -q {params.sample_name}); then
            bcftools view -s {params.sample_name} -Oz -o {output.imputation_vcf} {input.imputation_result}
        else
            touch {output.imputation_vcf}
        fi
    """

#vcf_dict = {}
#for panel in panels:
#    for id in seq_names:
#        vcf_dict[panel + id] = ["results/imputation/tmp/"+id+"/"+panel+"_chr"+str(chr)+".vcf.gz" for chr in chromosomes]

rule calculate_imputation_accuracy:
    input:
        imputation_vcf = expand("results/imputation/tmp/{id}/{panel}_chr{chr}.vcf.gz", chr = chromosome, allow_missing=True),
        chip_vcf = "results/chip/tmp/{id}/{id}.vcf.gz",
        afs = expand("data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt", chr = chromosome)
    output:
        r2 = "results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv"
    resources:
        mem_mb = 300000
    threads: 4
    run:
        chromosomes = chromosome
        vcfs = input.imputation_vcf
        mafs = input.afs
        chip = input.chip_vcf
        MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.95, 1])

        check_lst = vcfs + [chip]
        for i in check_lst: # If there is an empty file, report all imputation accuracy to be -1
            if os.path.getsize(i) == 0:
                r2 = -1*np.ones((2, np.size(MAF_ary) - 1))
                r2 = pd.DataFrame(r2.T, columns = ['Imputation Accuracy','Bin Count'], index = MAF_ary[1:])
                r2.index.name = 'MAF'
                r2.to_csv(output.r2, sep=',', mode='a')
                sys.exit(0)

        vcf = lcwgsus.multi_parse_vcf(chromosomes, vcfs)
        af = lcwgsus.multi_read_af(chromosomes, mafs)
        chip = lcwgsus.read_vcf(chip)
        chip = lcwgsus.drop_cols(chip, drop_lst = ['id', 'qual', 'filter','info','format'])
        r2 = lcwgsus.calculate_imputation_accuracy(vcf, chip, af)
        r2.to_csv(output.r2, sep=',', mode='a')

rule plot_imputation_accuracy:
    input:
        r2 = expand("results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv", id = seq_to_extract, panel = panels)
    output:
        graph = "graphs/imputation_vs_chip.png"
    resources:
        mem_mb = 30000
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

