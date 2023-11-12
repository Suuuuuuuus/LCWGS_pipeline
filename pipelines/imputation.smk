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
        vcf = f"results/imputation/vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem_mb = 30000
    params:
        threads = 1,
        nCores = 1
    resources:
        mem_mb = 30000
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        ## set a seed here, randomly, so can try to reproduce if it fails
        SEED=`echo $RANDOM`
        mkdir -p results/imputation/vcfs/regions/
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        nCores = {params.nCores}, \
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
        vcf = f"results/imputation/vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    output:
        vcf = f"results/imputation/vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz.output.RData"
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
        file="results/imputation/vcfs/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/imputation/vcfs/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule concat:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/imputation/vcfs/quilt.chr{{chr}}.vcf.gz"
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

        if [[ -d {params.rename_samples_file} && {params.rename_samples} == "True"]]
        then
            mv {output.vcf} tmp.{output.vcf}
            rm {output.vcf}.tbi
            bcftools reheader -o {output.vcf} -s {params.rename_samples_file} tmp.{output.vcf}
            tabix {output.vcf}
            rm tmp.{output.vcf}
        fi
    """

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
        if (bcftools query -l {input.chip_result} | grep -q {wildcards.id}); then
            bcftools view -s {params.chip_name} -Oz -o {output.chip_vcf} {input.chip_result}
        else
            touch {output.chip_vcf}
        fi
    """

rule get_imputation_vcf:
    input:
        imputation_result = "results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz"
    output:
        imputation_vcf = temp("results/imputation/tmp/{id}/{panel}_chr{chr}.vcf.gz"),
    params:
        seq_name = lambda wildcards: wildcards.id,
        sample_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Sample_Name'].values[0],
    resources:
        mem_mb = 30000
    shell: """
        if (bcftools query -l {input.imputation_result} | grep -q {wildcards.id}); then
            bcftools view -s {params.sample_name} -Oz -o {output.imputation_vcf} {input.imputation_>
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
        afs = expand("data/gnomAD_MAFs/gnomAD_MAF_afr_chr{chr}.txt", chr = chromosome)
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
        r2 = expand("results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv", id = seq_names, panel = panels)
    output:
        graph = "graphs/imputation_vs_chip.png"
    resources:
        mem_mb = 30000
    params:
        samples = seq_names,
        panels = panels
    run:
        samples = params.seq_names
        panels = params.panels
        dfs = []
        for i in panels:
            for j in samples:
                tmp = pd.read_csv("results/imputation/imputation_accuracy/"+j+"/"+i+"_imputation_accuracy.csv", sep = ',', dtype = {
                    'MAF': str,
                    'Imputation Accuracy': float,
                    'Bin Count': str
                }).iloc[2:,:]
                tmp['panel'] = i
                tmp['Bin Count'] = j
                tmp.columns = ['AF', 'corr', 'sample', 'panel']
                tmp['AF'] = tmp['AF'].shift(1).fillna('0') + '-' + tmp['AF']
                tmp['AF'] = tmp['AF'].astype("category")
                dfs.append(tmp)
        res = pd.concat(dfs).reset_index(drop = True)
        lcwgSus.plot_imputation_accuracy(res, plot_title = 'Imputation accuracy vs. ChIP, two reference panels', single_sample = False, save_fig = True, save_name = output.graph, outdir = None)