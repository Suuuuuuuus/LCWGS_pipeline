configfile: "pipelines/config.json"

import os
from os.path import exists
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("scripts")
import lcwgSus

# samples = pd.read_table(config['samples'], header = None, names = ['Code'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)
panels = config["panels"]

chip_to_extract = pd.read_table("chip.tsv", header = None).loc[:, 0].to_list()
seq_to_extract = sample_linker[sample_linker['Chip_Name'].isin(chip_to_extract)].Seq_Name.to_list()

chromosome = [i for i in range(1,23)]

# The followings are global parameters from `activate`:
QUILT_HOME = config["QUILT_HOME"]
ANALYSIS_DIR = config["ANALYSIS_DIR"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file) ## python is dumb

rule prepare_ref:
    input:
        json = "results/imputation/regions.json",
        hap = f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        recomb = f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem_mb = 30000
    params:
        threads = 8
    shell: """
        mkdir -p results/imputation/refs/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        nGen={NGEN}, \
        reference_haplotype_file="{input.hap}" ,\
        reference_legend_file="{input.legend}", \
        genetic_map_file="{input.recomb}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        output_file="{output.RData}")'
    """

rule quilt:
    input:
        bamlist = "results/imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem_mb = 30000
    resources:
        mem_mb = 30000
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    params:
        panel = PANEL_NAME
    shell: """
        ## set a seed here, randomly, so can try to reproduce if it fails
        SEED=`echo $RANDOM`
        mkdir -p "results/imputation/vcfs/{params.panel}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        bamlist="{input.bamlist}", \
        prepared_reference_filename="{input.RData}", \
        output_filename="{output.vcf}", \
        seed='${{SEED}}')'
    """
'''
rule quilt_info:
    input:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    output:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz.output.RData"
    params:
        threads = 1
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        R -f ${{QUILT_WRAP_HOME}}info.R --args {output.vcf}
    """
'''

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule concat:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
    resources:
        mem_mb = 30000
    params:
        threads = 1,
        input_string=get_input_vcfs_as_string,
        rename_samples = config["rename_samples"],
        rename_samples_file = config["rename_samples_file"]
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        if [ -e {output.vcf}.temp1.vcf.gz ]; then
            rm {output.vcf}.temp1.vcf.gz
        fi
        if [ -e {output.vcf}.temp2.vcf.gz ]; then
            rm {output.vcf}.temp2.vcf.gz
        fi

        bcftools concat \
        --ligate-force \
        --output-type z \
        --output {output.vcf}.temp1.vcf.gz \
        {params.input_string}

        gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS\t[%GT:%GP:%DS\t]\n' {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf
        bcftools sort -Oz -o {output.vcf} {output.vcf}.temp2.vcf
        tabix {output.vcf}
        rm {output.vcf}.temp*
    """
"""
        if [[ -d {params.rename_samples_file} && {params.rename_samples} == "True"]]
        then
            mv {output.vcf} tmp.{output.vcf}
            rm {output.vcf}.tbi
            bcftools reheader -o {output.vcf} -s {params.rename_samples_file} tmp.{output.vcf}
            tabix {output.vcf}
            rm tmp.{output.vcf}
        fi
"""

chip_to_extract = pd.read_table("chip.tsv", header = None).loc[:, 0].to_list()
seq_to_extract = sample_linker[sample_linker['Chip_Name'].isin(chip_to_extract)].Seq_Name.to_list()

rule get_chip_vcf:
    input:
        chip_result = "results/chip/filtered_snps.vcf.gz"
    output:
        chip_vcf = temp("results/chip/tmp/{id}/{id}.vcf.gz")
    params:
        seq_name = lambda wildcards: wildcards.id,
        chip_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Chip_Name'].values[0]
    resources:
        mem_mb = 30000
    shell: """
        if (bcftools query -l {input.chip_result} | grep -q {params.chip_name}); then
            bcftools view -s {params.chip_name} -Oz -o {output.chip_vcf} {input.chip_result}
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

        vcf = lcwgSus.multi_parse_vcf(chromosomes, vcfs)
        af = lcwgSus.multi_read_af(chromosomes, mafs)
        chip = lcwgSus.read_vcf(chip)
        chip = lcwgSus.drop_cols(chip, drop_lst = ['id', 'qual', 'filter','info','format'])
        r2 = lcwgSus.calculate_imputation_accuracy(vcf, chip, af)
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
        r2_fv = lcwgSus.read_r2(panels, samples_fv)
        r2_fv = lcwgSus.aggregate_r2(r2_fv)
        r2_mini = lcwgSus.read_r2(panels, samples_mini)
        r2_mini = lcwgSus.aggregate_r2(r2_mini)

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
