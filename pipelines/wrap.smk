configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

rule aggregate_results:
    input:
        samtools_dup_rate = rules.calculate_samtools_dup_rate.output.samtools_dup_rate_res,
	
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

rule extract_fastqc_dup_rate:
    input:
        zip1 = "results/fastqc/{id}_1_fastqc.zip",
        zip2 = "results/fastqc/{id}_2_fastqc.zip"
    output:
        txt = temp("results/fastqc/tmp/{id}.txt")
    params:
        tmpdir = "results/fastqc/"
    shell: """
        unzip "{params.tmpdir}{wildcards.id}_1_fastqc.zip" -d {params.tmpdir}
        unzip "{params.tmpdir}{wildcards.id}_2_fastqc.zip" -d {params.tmpdir}
        dedup1=$(grep "Total Deduplicated Percentage" "{params.tmpdir}{wildcards.id}_1_fastqc/fastqc_data.txt" | cut -f 2)
        dedup2=$(grep "Total Deduplicated Percentage" "{params.tmpdir}{wildcards.id}_2_fastqc/fastqc_data.txt" | cut -f 2)
        echo "{wildcards.id} $dedup1 $dedup2" > {output.txt}
        rm -r {params.tmpdir}*{wildcards.id}*fastqc
    """

rule aggregate_fastqc_dup_rate:
    input:
        files = expand("results/fastqc/tmp/{id}.txt", id = ids_1x_all)
    output:
        fastqc_dup_rate = temp("results/fastqc/duplication_rate_fastqc_tmp.txt")
    params:
        tmpdir = "results/fastqc/tmp/"
    shell: """
        cat {input.files} >> {output.fastqc_dup_rate}
        rm -r {params.tmpdir}
    """

rule calculate_fastqc_dup_rate:
    input:
        fastqc_dup_rate = rules.aggregate_fastqc_dup_rate.output.fastqc_dup_rate
    output:
        fastqc_dup_rate_res = "results/fastqc/duplication_rate_fastqc.txt"
    shell: """
        while IFS=$' ' read -r col1 col2 col3; do
            result=$(echo "scale=4; (1-($col3+$col2)/200)" | bc)
            formatted_result=$(printf "%06.4f" $result)
            echo -e "$col1\t$formatted_result" >> {output.fastqc_dup_rate_res}
        done < {input.fastqc_dup_rate}
    """
