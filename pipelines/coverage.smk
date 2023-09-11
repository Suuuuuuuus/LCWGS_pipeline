configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
chromosome = [i for i in range(1,23)]

subsample_coverage = config['subsample_depth']

rule compute_bedgraph:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph.txt"
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
        ss_bedgraph = "results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.txt"
    resources: mem_mb = 50000
    shell: """
        bedtools genomecov -ibam {input.ss_bam} -bga | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.ss_bedgraph}
    """

rule calculate_ss_cumsum_coverage:
    input:
        ss_bedgraphs = "results/coverage/subsampled_bedgraphs/{id}_subsampled_bedgraph.txt",
        code = "scripts/calculate_ss_cumsum_coverage.py"
    output:
        cumsum_ary = "results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt"
    params:
        num_coverage = 10, # Specify the length of the x-axis
        avg_coverage = subsample_coverage # Specify the Poisson expectation loc parameter
    resources: mem_mb = 5000
    shell: """
        python {input.code} {input.ss_bedgraphs} {params.num_coverage} {params.avg_coverage}
    """

rule plot_subsample_coverage:
    input:
        code = "scripts/plot_subsample_coverage.py",
        cumsum_ary = expand("results/coverage/subsampled_bedgraphs/{id}_cumsum_ary.txt", id = ids_1x_all)
    output:
        graph = "graphs/fig8_prop_genome_at_least_coverage.png"
    resources: mem_mb = 5000
    params:
        num_coverage = 10, # Specify the length of the x-axis
        avg_coverage = subsample_coverage # Specify the Poisson expectation loc parameter
    shell: """
        python {input.code} {params.num_coverage} {params.avg_coverage}
    """
'''
rule calculate_per_bin_coverage_1x:
    input:
        script = "scripts/calculate_per_bin_coverage.py",
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph.txt"
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
'''

rule samtools_coverage:
    input:
        bam = "data/bams/{id}.bam"
    output:
        per_chromosome_coverage = "results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt"
    shell: """
        samtools coverage {input.bam} | sed -n '2,23p' > {output.per_chromosome_coverage}
    """

rule samtools_ss_coverage:
    input:
        ss_bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        per_chromosome_ss_coverage = "results/coverage/per_chromosome_ss_coverage/{id}_per_chromosome_ss_coverage.txt"
    shell: """
        samtools coverage {input.ss_bam} | sed -n '2,23p' > {output.per_chromosome_ss_coverage}
    """

rule calculate_ss_uncoverage_rate:
    input:
        per_chromosome_ss_coverage = rules.samtools_ss_coverage.output.per_chromosome_ss_coverage
    output:
        ss_uncoverage_rate = temp("results/coverage/per_chromosome_ss_coverage/{id}_ss_uncoverage_rate.txt"),
    shell: """
        total=$(cut -f3 {input.per_chromosome_ss_coverage} | paste -sd+ | bc)
        covered=$(cut -f5 {input.per_chromosome_ss_coverage} | paste -sd+ | bc)
        result=$(echo "scale=4; (1-$covered/$total)" | bc)
        echo "{wildcards.id}\t$result" > {output.ss_uncoverage_rate}
    """

rule aggregate_ss_uncoverage_rate:
    input:
        files = expand("results/coverage/per_chromosome_ss_coverage/{id}_ss_uncoverage_rate.txt", id = ids_1x_all)
    output:
        ss_uncoverage_rate = "results/coverage/per_chromosome_ss_coverage/ss_uncoverage_rate.txt"
    shell: """
        cat {input.files} >> {output.ss_uncoverage_rate}
    """

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
rule plot_uncoverage_rate:
    input:
        code = "scripts/plot_uncoverage_rate.py",
	uncoverage_rate = rules.aggregate_uncoverage_rate.output.uncoverage_rate
    output:
        graph = "graphs/uncoverage_rate.png"
    resources: mem_mb = 50000
    shell: """
	python {input.code}
    """

rule plot_per_bin_coverage:
    input:
        code = "scripts/plot_per_bin_coverage.py",
	bases = expand("results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt", id = ids_1x_all, allow_missing=True)
    output:
        graph = "graphs/fig6_per_bin_coverage_chr{chr}.png"
    params:
        repeat_mask_bin_size = 100000,
        ylim = 20
    shell: """
        python {input.code} {wildcards.chr} {params.repeat_mask_bin_size} {params.ylim}
    """
'''
