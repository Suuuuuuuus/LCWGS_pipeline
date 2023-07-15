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

rule prepare_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = ids_1x_all)
    output:
        bamlist = "results/imputation/bamlist.txt"
    params:
        threads=1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        mkdir -p {ANALYSIS_DIR}
        ls /well/band/users/rbx225/GGVP/data/bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f ${{QUILT_HOME}}scripts/make_b38_recomb_map.R --args "./" {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/imputation_refs/ggvp.chr{{chr}}.vcf.gz",
        tbi = f"data/imputation_refs/ggvp.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp("results/imputation/refs/tmp.ggvp.chr{chr}.vcf.gz"),
        hap = temp("results/imputation/refs/ggvp.chr{chr}.hap.gz"),
        legend = temp("results/imputation/refs/ggvp.chr{chr}.legend.gz"),
        samples = temp("results/imputation/refs/ggvp.chr{chr}.samples")
    params:
        threads=1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        mkdir -p results/imputation/refs/
        bcftools view --output-file {output.tmp_vcf} --output-type z --min-alleles 2 --max-alleles 2 --types snps {input.vcf}
        tabix {output.tmp_vcf}
        bcftools convert --haplegendsample results/imputation/refs/ggvp.chr{wildcards.chr} {output.tmp_vcf}
        rm {output.tmp_vcf} {output.tmp_vcf}.tbi
    """

rule determine_chunks:
    input:
        legend = expand("results/imputation/refs/ggvp.chr{chr}.legend.gz", chr = chromosome)
    output:
        json = "results/imputation/regions.json"
    shell: """
        R -f scripts/determine_chunks.R
    """