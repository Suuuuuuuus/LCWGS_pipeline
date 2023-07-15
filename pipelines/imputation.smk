configfile: "pipelines/config.json"

from os.path import exists
import json
import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
chromosome = [i for i in range(1,23)]

# The followings are global parameters from `activate`:
QUILT_HOME = "/well/band/users/rbx225/software/QUILT/"
ANALYSIS_DIR = "/well/band/users/rbx225/GGVP/results/imputation/"
RECOMB_POP="ACB"
NGEN=100
WINDOWSIZE=5000000
BUFFER=1000000

file="results/imputation/regions.json"
with open(file) as json_file:
    REGIONS = json.load(json_file) ## python is dumb

rule prepare_ref:
    input:
        json = "results/imputation/regions.json",
        hap = "results/imputation/refs/ggvp.chr{chr}.hap.gz",
        legend = "results/imputation/refs/ggvp.chr{chr}.legend.gz",
        recomb = f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    params:
        threads = 1
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

rule imputation_accuracy_NA12878:
        input:
                script = "scripts/calculate_imputation_accuracy_NA12878_v2.py"
        output:
                graph = "graphs/NA12878_imputation_accuracy.png"
        resources:
                mem_mb = 300000
        threads:
                16
        shell: """
                python {input.script}
        """

rule imputation_accuracy_lcwgs:
        input:
                script = "scripts/calculate_imputation_accuracy_lcwgs_v2.py"
        output:
                graph = "graphs/lcwgs_imputation_accuracy.png"
        resources:
                mem_mb = 300000
        threads:
                16
        shell: """
                python {input.script}
        """

'''
rule extract_gnomAD_MAF:
	input:	
		gnomAD_vcf = "data/gnomAD_vcf/gnomad.genomes.v3.1.2.sites.chr{chr}.vcf.bgz"
	output:
		gnomAD_MAF = "results/variant_calling/gnomAD_MAFs/gnomAD_MAF_chr{chr}.txt"
	threads: 4
	shell: """
		/bin/bash -c "zgrep -v '#' {input.gnomAD_vcf} | \
		awk -v "OFS=\t" '{$8=$8;sub(/^.*AF=/, "", $8); sub(/;.*/, "", $8); print $1,$2,$4,$5,$8}' | \
		awk '$5 !~ /^[[:alpha:]]/' | \
		awk '!(length($3)>1 || length($4)>1)' | \
		awk '!($5==0.00000)' > {output.gnomAD_MAF}"
	"""
'''
