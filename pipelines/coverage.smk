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

rule compute_bedgraph:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph.bed"
    resources: mem_mb = 50000
    shell: """
        bedtools genomecov -ibam {input.bam} -bga | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
    """

rule compute_subsampled_bedgraph:
    input:
        ss_bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        ss_bedgraph = "results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.bed"
    resources: mem_mb = 50000
    shell: """
        bedtools genomecov -ibam {input.ss_bam} -bga | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.ss_bedgraph}
    """

rule calculate_ss_cumsum_coverage:
    input:
        ss_bedgraph = "results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.bed"
    output:
        cumsum_ary = "results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt"
    params:
        num_coverage = config["num_coverage"], 
        rm_bed_regions = config["rm_bed_regions"],
        bed_regions = config["bed_regions"]
    resources:
        mem_mb = 20000
    threads: 11
    run:
        chromosomes = [i for i in range(1,23)]
        cols = ['chr', 'start', 'end', 'cov']
        cov = pd.read_csv(input.ss_bedgraph, sep = '\t', header = None, names = cols)
        cov['chr'] = cov['chr'].str.extract(r'(\d+)').astype(int)
        cov = cov.sort_values(by = ['chr', 'start'])
        if params.rm_bed_regions:
            region = pd.read_csv(params.bed_regions, sep = '\t', header = None, names = cols[:-1])
            region = region[region['chr'].isin(['chr' + str(i) for i in range(1,23)])]
            region['chr'] = region['chr'].str.extract(r'(\d+)').astype(int)
            region = region.sort_values(by = ['chr', 'start'])
            bed = lcwgSus.multi_subtract_bed(chromosomes, lcwgSus.file_to_list(cov), lcwgSus.file_to_list(region))
        else:
            bed = cov
        cumsum_ary = lcwgSus.calculate_ss_cumsum_coverage(bed, num_coverage = params.num_coverage)
        np.savetxt(output.cumsum_ary, cumsum_ary, newline = ' ', fmt = '%1.7f')

rule plot_sequencing_skew:
    input:
        cumsum_ary = expand("results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt", id = ids_1x_all)
    output:
        graph_seq_skew = "graphs/prop_genome_at_least_coverage.png"
    resources: mem_mb = 5000
    params:
        num_coverage = config["num_coverage"], # Specify the length of the x-axis
        avg_coverage = config["subsample_depth"] # Specify the Poisson expectation loc parameter
    run:
        ary_lst = []
        for i in input.cumsum_ary:
            ary_lst.append(np.loadtxt(i))
        lcwgSus.plot_sequencing_skew(ary_lst, params.avg_coverage, num_coverage = params.num_coverage, save_fig = True)

'''
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
'''
rule calculate_uncoverage_rate:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = temp("results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"),
        uncoverage_rate = temp("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt")
    params:
        access_bed = config['access_bed']
    resources:
        mem_mb = 30000
    shell: """
        bedtools genomecov -ibam {input.bam} -bg | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
        result=$(bedtools coverage -a {params.access_bed} -b {output.bedgraph} -hist | grep all | head -n 1 | cut -f5)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
    """

rule aggregate_uncoverage_rate:
    input:
        files = expand("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt", id = ids_1x_all)
    output:
        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt"
    shell: """
        cat {input.files} >> {output.uncoverage_rate}
    """

rule calculate_ss_uncoverage_rate:
    input:
        ss_bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        bedgraph = temp("results/coverage/bedgraphs/{id}_ss_bedgraph_nozero.bed"),
        uncoverage_rate = temp("results/coverage/per_chromosome_ss_coverage/{id}_ss_uncoverage_rate.txt")
    params:
        access_bed = config['access_bed']
    resources:
        mem_mb = 30000
    shell: """
        bedtools genomecov -ibam {input.ss_bam} -bg | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
        result=$(bedtools coverage -a {params.access_bed} -b {output.bedgraph} -hist | grep all | head -n 1 | cut -f5)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
    """

rule aggregate_ss_uncoverage_rate:
    input:
        files = expand("results/coverage/per_chromosome_ss_coverage/{id}_ss_uncoverage_rate.txt", id = ids_1x_all)
    output:
        ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt"
    shell: """
        cat {input.files} >> {output.ss_uncoverage_rate}
    """
'''
rule calculate_avg_coverage:
    input:
        per_chromosome_coverage = rules.samtools_coverage.output.per_chromosome_coverage
    output:
        avg_coverage = temp("results/coverage/tmp/{id}_avg_coverage.txt")
    shell: """
        total=$(cut -f3 {input.per_chromosome_coverage} | paste -sd+ | bc)
        sum_product=$(awk '{{ sum += $3 * $7 }} END {{ printf "%.2f", sum }}' {input.per_chromosome_coverage})
        result=$(echo "scale=4; ($sum_product/$total)" | bc)
        echo "{wildcards.id}\t$result" > {output.avg_coverage}
    """

rule aggregate_avg_coverage:
    input:
        files = expand("results/coverage/tmp/{id}_avg_coverage.txt", id = ids_1x_all)
    output:
        avg_coverage = "results/coverage/per_sample_coverage.txt"
    shell: """
        cat {input.files} >> {output.avg_coverage}
    """
'''
