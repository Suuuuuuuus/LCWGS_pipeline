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
        outputdir = "results/fastqc/"
    shell: """
        mkdir -p results
        mkdir -p results/fastqc
        mkdir -p {params.outputdir}
        fastqc -q -o {params.outputdir} {input.fastq1} {input.fastq2}
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
        outputdir = "results/fastqc/"
    shell: """
        mkdir -p results
        mkdir -p results/fastqc
        mkdir -p {params.outputdir}
        fastqc -q -o {params.outputdir} {input.fastq1} {input.fastq2}
    """

rule multiqc_lc:
    input:
        expand("results/fastqc/{id}_1_fastqc.html",id = ids_1x_all),
        expand("results/fastqc/{id}_2_fastqc.html",id = ids_1x_all),
        expand("results/fastqc/{id}_1_fastqc.zip",id = ids_1x_all),
        expand("results/fastqc/{id}_2_fastqc.zip",id = ids_1x_all)
    output:
        "results/fastqc/multiqc_report_lc.html",
        directory("results/fastqc/multiqc_data")
    threads: 1
    resources:
        mem = '10G'
    params:
        "results/fastqc/"
    shell:
        "multiqc {params} --interactive -o {params}"

rule multiqc_hc:
    input:
        expand("results/fastqc/{id}_1_fastqc.html",id = samples_hc),
        expand("results/fastqc/{id}_2_fastqc.html",id = samples_hc),
        expand("results/fastqc/{id}_1_fastqc.zip",id = samples_hc),
        expand("results/fastqc/{id}_2_fastqc.zip",id = samples_hc)
    output:
        "results/fastqc/multiqc_report.html",
        directory("results/fastqc/multiqc_data")
    threads: 1
    resources:
        mem = '10G'
    params:
        "results/fastqc/"
    shell:
        "multiqc {params} --interactive -o {params}"
