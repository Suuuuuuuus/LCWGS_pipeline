configfile: "pipelines/config.json"

import pandas as pd
import numpy as np
import sys
sys.path.append("scripts")
import lcwgSus

config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
chromosome = [i for i in range(1,23)]

subsample_coverage = config['subsample_depth']
rm_bed_regions = config['rm_bed_regions']
bed_regions = config['bed_regions']

rule calculate_per_bin_coverage_1x:
    input:
        script = "scripts/calculate_per_bin_coverage.py",
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph.bed"
    output:
        coordinates = "results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt",
        bases = "results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt"
    params:
        indir = "results/coverage/bedgraphs/",
        outdir = "results/coverage/per_bin_coverage/1x/",
        bin_size = 10000
    threads: 8
    resources: mem_mb = 50000
    shell: """
        python {input.script} {wildcards.id} {wildcards.chr} {params.indir} {params.outdir} {params.bin_size}
    """

rule calculate_per_bin_coverage_20x:
    input:
        script = "scripts/calculate_per_bin_coverage.py"
    output:
        coordinates = "results/coverage/per_bin_coverage/20x/{id_20x}_chr{chr}_coordinate.txt",
        bases = "results/coverage/per_bin_coverage/20x/{id_20x}_chr{chr}_base.txt"
    params:
        indir = "results/coverage/bedgraphs/",
        outdir = "results/coverage/per_bin_coverage/20x/",
        bin_size = 10000
    # threads: 16
    # resources: mem_mb = 40000
    shell: """
        python {input.script} {wildcards.id_20x} {wildcards.chr} {params.indir} {params.outdir} {params.bin_size}
    """

rule samtools_coverage:
    input:
        bam = "data/bams/{id}.bam"
    output:
        per_chromosome_coverage = "results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt"
    shell: """
        samtools coverage {input.bam} | sed -n '2,23p' > {output.per_chromosome_coverage}
    """

rule calculate_uncoverage_rate:
    input:
        per_chromosome_coverage = rules.samtools_coverage.output.per_chromosome_coverage
    output:
        uncoverage_rate = temp("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt")
    shell: """
	total=$(cut -f3 {input.per_chromosome_coverage} | paste -sd+ | bc)
        covered=$(cut -f5 {input.per_chromosome_coverage} | paste -sd+ | bc)
        result=$(echo "scale=4; (1-$covered/$total)" | bc)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
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