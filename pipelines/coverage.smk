configfile: "pipelines/config.json"

import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

subsample_coverage = config['subsample_depth']
rm_bed_regions = config['rm_bed_regions']
bed_regions = config['bed_regions']

rule compute_bedgraph_nozero:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    resources: mem_mb = 50000
    shell: """
        bedtools genomecov -ibam {input.bam} -bg | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
    """

rule calculate_avg_coverage:
    input:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    output:
        access_coverage = temp("results/coverage/tmp/{id}_access_coverage.txt"),
        avg_coverage = temp("results/coverage/tmp/{id}_avg_coverage.txt")
    params:
        access_bed = config['access_bed'],
        access_bed_length = 2526390487
    shell: """
        bedtools intersect -a {input.bedgraph} -b {params.access_bed} -wb | cut -f1-4 > {output.access_coverage}
        #total=$
        sum_product=$(awk '{{ sum += ($3-$2)*$4 }} END {{ print sum }}' {output.access_coverage})
        result=$(echo "scale=4; ($sum_product/{params.access_bed_length})" | bc)
        echo "{wildcards.id}\t$result" > {output.avg_coverage}
    """

rule aggregate_avg_coverage:
    input:
        files = expand("results/coverage/tmp/{id}_avg_coverage.txt", id = samples_lc)
    output:
        avg_coverage = "results/coverage/per_sample_coverage.txt"
    shell: """
        cat {input.files} >> {output.avg_coverage}
    """

rule aggregate_hc_avg_coverage:
    input:
        files = expand("results/coverage/tmp/{id}_avg_coverage.txt", id = samples_hc)
    output:
        hc_coverage = "results/coverage/hc_coverage.txt"
    shell: """
        cat {input.files} >> {output.hc_coverage}
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
            bed = lcwgsus.multi_subtract_bed(chromosomes, lcwgsus.file_to_list(cov), lcwgsus.file_to_list(region))
        else:
            bed = cov
        cumsum_ary = lcwgsus.calculate_ss_cumsum_coverage(bed, num_coverage = params.num_coverage)
        np.savetxt(output.cumsum_ary, cumsum_ary, newline = ' ', fmt = '%1.7f')

rule plot_sequencing_skew:
    input:
        cumsum_ary = expand("results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt", id = samples_lc)
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
        lcwgsus.plot_sequencing_skew(ary_lst, params.avg_coverage, num_coverage = params.num_coverage, save_fig = True)


rule calculate_uncoverage_rate:
    input:
        bam = "data/bams/{id}.bam",
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    output:
        uncoverage_rate = temp("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt")
    params:
        access_bed = config['access_bed']
    resources:
        mem_mb = 30000
    shell: """
        result=$(bedtools coverage -a {params.access_bed} -b {input.bedgraph} -hist | grep all | head -n 1 | cut -f5)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
    """

rule aggregate_uncoverage_rate:
    input:
        files = expand("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt", id = samples_lc)
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
        files = expand("results/coverage/per_chromosome_ss_coverage/{id}_ss_uncoverage_rate.txt", id = samples_lc)
    output:
        ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt"
    shell: """
        cat {input.files} >> {output.ss_uncoverage_rate}
    """