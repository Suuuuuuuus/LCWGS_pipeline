configfile: "pipelines/config.json"

from os.path import exists
import json
import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
# chromosome = [i for i in range(1,23)]
chromosome = [11]

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
    resources:
        mem_mb = 30000
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

rule quilt:
    input:
        bamlist = "results/imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/imputation/vcfs/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem_mb = 30000
    params:
        threads = 1
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
        input_string=get_input_vcfs_as_string
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        bcftools concat \
        --ligate \
        --output-type z \
        --output {output.vcf}.temp1.vcf.gz \
        {params.input_string}

        gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS:PS\t[%GT:%GP:%DS:%PS\t]\n' {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf
        bgzip {output.vcf}.temp2.vcf
        tabix {output.vcf}.temp2.vcf.gz

        mv {output.vcf}.temp2.vcf.gz {output.vcf}
        mv {output.vcf}.temp2.vcf.gz.tbi {output.vcf}.tbi
        rm {output.vcf}.temp1.vcf.gz
    """

'''
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
