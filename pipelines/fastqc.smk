configfile: "pipelines/config.json"

import pandas as pd
ids_1x_all = list(pd.read_table(config["samples"], header = None, names = ['Code'])['Code'].values)

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
        mkdir -p results
        mkdir -p results/fastqc
        mkdir -p {params.outdir}
        fastqc -q -o {params.outdir} {input.fastq1} {input.fastq2}
    """

rule multiqc_lc:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html",id = ids_1x_all),
        html2 = expand("results/fastqc/{id}_2_fastqc.html",id = ids_1x_all),
        zip = expand("results/fastqc/{id}_{read}_fastqc.zip",id = ids_1x_all, read = ['1', '2'])
    output:
        html = "results/fastqc/multiqc_lc/multiqc_report.html",
        directory("results/fastqc/multiqc_lc")
    threads: 1
    resources:
        mem = '10G'
    params:
        outdir = "results/fastqc/multiqc_lc"
    shell:
        "multiqc {input.zip} --interactive -o {params.outdir}"

rule multiqc_lc:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html",id = samples_hc),
        html2 = expand("results/fastqc/{id}_2_fastqc.html",id = samples_hc),
        zip = expand("results/fastqc/{id}_{read}_fastqc.zip",id = samples_hc, read = ['1', '2'])
    output:
        html = "results/fastqc/multiqc_hc/multiqc_report.html",
        directory("results/fastqc/multiqc_hc")
    threads: 1
    resources:
        mem = '10G'
    params:
        outdir = "results/fastqc/multiqc_hc"
    shell:
        "multiqc {input.zip} --interactive -o {params.outdir}"
