configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

rule extract_samtools_dup_rate:
    input:
        bam = "data/bams/{id}.bam"
    output:
        txt = temp("results/dup_rate/tmp/{id}.txt")
    params:
        tmpdir = "results/dup_rate/tmp/"
    shell: """
        samtools sort -n -o "{params.tmpdir}{wildcards.id}_sorted.bam" -O BAM {input.bam}
        samtools fixmate -m "{params.tmpdir}{wildcards.id}_sorted.bam" "{params.tmpdir}{wildcards.id}_fixmate.bam"
        samtools sort -o "{params.tmpdir}{wildcards.id}_position_sorted.bam" "{params.tmpdir}{wildcards.id}_fixmate.bam"
        samtools markdup "{params.tmpdir}{wildcards.id}_position_sorted.bam" "{params.tmpdir}{wildcards.id}_markdup.bam"
        samtools flagstat "{params.tmpdir}{wildcards.id}_markdup.bam" > "{params.tmpdir}{wildcards.id}_result.txt"

        total=$(cat "{params.tmpdir}{wildcards.id}_result.txt" | cut -f1 -d ' '| head -n 1)
        dup=$(cat "{params.tmpdir}{wildcards.id}_result.txt" | cut -f1 -d ' '| head -n 5 | tail -n 1)
        echo "{wildcards.id} $total $dup" > "{output.txt}"

        rm {params.tmpdir}{wildcards.id}*.bam
    """

rule aggregate_samtools_dup_rate:
    input:
        files = expand("results/dup_rate/tmp/{id}.txt", id = ids_1x_all)
    output:
        samtools_dup_rate = temp("results/dup_rate/duplication_rate_samtools_tmp.txt")
    params:
        tmpdir = "results/dup_rate/tmp/"
    shell: """
        cat {input.files} >> {output.samtools_dup_rate}
        rm -r {params.tmpdir}
    """

rule calculate_samtools_dup_rate:
    input:
        samtools_dup_rate = rules.aggregate_samtools_dup_rate.output.samtools_dup_rate
    output:
        samtools_dup_rate_res = "results/dup_rate/duplication_rate_samtools.txt"
    shell: """
        while IFS=$' ' read -r col1 col2 col3; do
            result=$(echo "scale=4; $col3/$col2" | bc)
            formatted_result=$(printf "%06.4f" $result)
            echo -e "$col1\t$formatted_result" >> {output.samtools_dup_rate_res}
        done < {input.samtools_dup_rate}
    """

rule plot_samtools_dup_rate:
    input:
        script = "scripts/plot_samtools_duplication_rate.py",
	samtools_dup_rate = rules.calculate_samtools_dup_rate.output.samtools_dup_rate_res
    output:
        graph_samtools_dup_rate = "graphs/samtools_duplication_rate.png"
    shell: """
        python {input.script}
    """
