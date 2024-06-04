configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])

if config['clean_fastq']:
    ruleorder: fastqc_alt > fastqc
else:
    ruleorder: fastqc > fastqc_alt

rule fastqc:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        html1 = "results/fastqc/{id}_1_fastqc.html",
        html2 = "results/fastqc/{id}_2_fastqc.html",
        zip1 = "results/fastqc/{id}_1_fastqc.zip",
        zip2 = "results/fastqc/{id}_2_fastqc.zip"
    params:
        outdir = "results/fastqc/"
    shell: """
        mkdir -p results
        mkdir -p results/fastqc
        mkdir -p {params.outdir}
        fastqc -q -o {params.outdir} {input.fastq1} {input.fastq2}
    """

rule fastqc_alt:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz"
    output:
        html1 = "results/fastqc/{id}_1_fastqc.html",
        html2 = "results/fastqc/{id}_2_fastqc.html",
        zip1 = "results/fastqc/{id}_1_fastqc.zip",
        zip2 = "results/fastqc/{id}_2_fastqc.zip"
    params:
        outdir = "results/fastqc/"
    shell: """
        mkdir -p {params.outdir}
        fastqc -q -o {params.outdir} {input.fastq1} {input.fastq2}
    """

rule multiqc_lc:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html",id = samples_lc),
        html2 = expand("results/fastqc/{id}_2_fastqc.html",id = samples_lc),
        zip = expand("results/fastqc/{id}_{read}_fastqc.zip",id = samples_lc, read = ['1', '2'])
    output:
        html = "results/fastqc/multiqc_lc/multiqc_report.html"
    threads: 1
    resources:
        mem = '10G'
    params:
        outdir = "results/fastqc/multiqc_lc"
    shell:
        "multiqc {input.zip} --interactive -o {params.outdir}"

samples_fv = read_tsv_as_lst("data/sample_tsvs/fv_gm_names.tsv")

rule get_per_base_error_rate:
    input:
        zip = "results/fastqc/{id}_{read}_fastqc.zip"
    output:
        tsv = temp("results/fastqc/{id}_{read}_per_base_error_rate.tsv")
    threads: 1
    resources:
        mem = '10G'
    params:
        outdir = "results/fastqc/tmp/{id}_read{read}/",
        data = "{id}_{read}_fastqc/fastqc_data.txt"
    shell: """
        mkdir -p {params.outdir}

        unzip -j {input.zip} {params.data} -d {params.outdir}
        sed -n '15,52p' {params.outdir}fastqc_data.txt | cut -f2 > {output.tsv}
    """

def convert_phred_to_prob(num):
    return 10**(-0.1*float(num))

rule calculate_average_per_base_error_rate:
    input:
        tsv1 = expand("results/fastqc/{id}_1_per_base_error_rate.tsv", id = samples_fv),
        tsv2 = expand("results/fastqc/{id}_2_per_base_error_rate.tsv", id = samples_fv)
    output:
        tsv = "results/fastqc/per_base_error_rate.tsv"
    threads: 1
    resources:
        mem = '10G'
    run:
        read1 = []
        read2 = []

        for i in input.tsv1:
            ary1 = read_tsv_as_lst(i)
            read1.append([convert_phred_to_prob(ary1[0]), convert_phred_to_prob(ary1[-1])])
        for j in input.tsv2:
            ary2 = read_tsv_as_lst(j)
            read2.append([convert_phred_to_prob(ary2[0]), convert_phred_to_prob(ary2[-1])])

        read1 = np.array(read1)
        read2 = np.array(read2)
        r1_b1, r1_b151 = read1.mean(axis = 0)
        r2_b1, r2_b151 = read2.mean(axis = 0)

        res = pd.DataFrame({'start': [r1_b1, r2_b1], 'end': [r1_b151, r2_b151]})
        print(res)
        res.to_csv(output.tsv, sep = '\t', index = False, header = True)
