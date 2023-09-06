configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

rule extract_samtools_dup_rate:
    input:
        bam = "data/bams/{id}.bam"
    output:
        txt = temp("results/dup_rate/tmp/{id}.txt")
    resources: mem_mb = 50000
    params:
        tmpdir = "results/dup_rate/tmp/"
    shell: """
        samtools sort -n -o "{params.tmpdir}{wildcards.id}_sorted.bam" -O BAM {input.bam}
        samtools fixmate -m "{params.tmpdir}{wildcards.id}_sorted.bam" "{params.tmpdir}{wildcards.id}_fixmate.bam"
        rm "{params.tmpdir}{wildcards.id}_sorted.bam"
        samtools sort -o "{params.tmpdir}{wildcards.id}_position_sorted.bam" "{params.tmpdir}{wildcards.id}_fixmate.bam"
        rm "{params.tmpdir}{wildcards.id}_fixmate.bam"
        samtools markdup "{params.tmpdir}{wildcards.id}_position_sorted.bam" "{params.tmpdir}{wildcards.id}_markdup.bam"
        rm "{params.tmpdir}{wildcards.id}_position_sorted.bam"
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

rule calculate_fragment_size:
    input:
        bam = "data/bams/{id}.bam"
    output:
        txt = temp("results/fragment_size/{id}/{id}.txt")
    shell: """
        mkdir -p results/fragment_size/
        samtools stat {input.bam} | grep ^IS | cut -f 2,3 > results/fragment_size/{wildcards.id}/tmp.txt
        total=$(cut -f2 results/fragment_size/{wildcards.id}/tmp.txt | paste -sd+ | bc)
        sum_product=$(awk '{{ sum += $1 * $2 }} END {{ printf "%.2f", sum }}' results/fragment_size/{wildcards.id}/tmp.txt)
        result=$(echo "scale=4; ($sum_product/$total)" | bc)
        echo "{wildcards.id}\t$result" > {output.txt}
    """

rule aggregate_fragment_size:
    input:
        files = expand("results/fragment_size/{id}/{id}.txt", id = ids_1x_all)
    output:
        fragment_size = "results/fragment_size/fragment_size.txt"
    shell: """
        cat {input.files} >> {output.fragment_size}
    """

rule calculate_proportion_ss_fragment_size:
    input:
        ss_bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        fragment_size = temp("results/fragment_size/{id}/fragment_size.txt"),
        txt = temp("results/fragment_size/{id}/{id}_proportion.txt")
    params:
        threshold = 500
    shell: """
        samtools view -f 66 -F 256 {input.ss_bam} | cut -f9 > {output.fragment_size}
        #file_path={output.fragment_size}
        propn_count=$(awk '{{ if ($1 > {params.threshold} || $1 < -{params.threshold}) count++ }} END {{ print count }}' {output.fragment_size})
        total_count=$(wc -l < {output.fragment_size})
        proportion=$(awk "BEGIN {{ print $propn_count / $total_count }}")
        echo "{wildcards.id}\t$proportion" > {output.txt}
    """

rule aggregate_proportion_ss_fragment_size:
    input:
        files = expand("results/fragment_size/{id}/{id}_proportion.txt", id = ids_1x_all)
    output:
        proportion_ss_fragment_size = "results/fragment_size/proportion_ss_fragment_size.txt"
    shell: """
        cat {input.files} >> {output.proportion_ss_fragment_size}
    """

rule calculate_proportion_fragment_size:
    input:
        bam = "data/bams/{id}.bam"
    output:
        fragment_size = temp("results/fragment_size/{id}/fragment_size_whole.txt"),
        txt = temp("results/fragment_size/{id}/{id}_proportion_whole.txt")
    params:
        threshold = 302
    shell: """
	samtools view -f 66 -F 256 {input.bam} | cut -f9 > {output.fragment_size}
        file_path={output.fragment_size}
        propn_count=$(awk '{{ if ($1 < {params.threshold} && $1 > -{params.threshold}) count++ }} END {{ print count }}' "$file_path")
        total_count=$(wc -l < "$file_path")
        proportion=$(awk "BEGIN {{ print $propn_count / $total_count }}")
	echo "{wildcards.id}\t$proportion" > {output.txt}
    """

rule aggregate_proportion_fragment_size:
    input:
        files = expand("results/fragment_size/{id}/{id}_proportion_whole.txt", id = ids_1x_all)
    output:
        proportion_fragment_size = "results/fragment_size/proportion_fragment_size.txt"
    shell: """
	cat {input.files} >> {output.proportion_fragment_size}
    """

rule calculate_fragment_overlap:
    input:
        bam = "data/bams/{id}.bam"
    output:
        fragment_overlap = temp("results/fragment_size/{id}/fragment_overlap.txt")
    params:
        threshold = 302
    shell: """
        total_overlap=$(samtools view -f 66 -F 256 {input.bam} | cut -f9 | \
        awk '{{$1 = ($1 < 0 ? -$1 : $1); $1 = ({params.threshold} - $1); $1 = (0 > $1 ? 0 : $1); print}}' | \
        paste -sd+ | bc)
        echo "{wildcards.id}\t$total_overlap" > {output.fragment_overlap}
    """

rule aggregate_fragment_overlap:
    input:
        files = expand("results/fragment_size/{id}/fragment_overlap.txt", id = ids_1x_all)
    output:
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"
    shell: """
        cat {input.files} >> {output.fragment_overlap}
    """
