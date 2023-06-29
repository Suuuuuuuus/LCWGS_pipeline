configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

rule classify_kmers:
    input:
        ss_fastq = "data/subsampled_fastq/{id}_subsampled_{read}.fastq",
        jf = "data/jellyfish/NA12878_31mer.jf"
    output:
        jf_read = "results/kmer/{id}/read{read}/{id}_read{read}.tsv",
        jf_quality = "results/kmer/{id}/read{read}/{id}_quality{read}.tsv",
        jf_position = "results/kmer/{id}/read{read}/{id}_position{read}.tsv"
    params:
        read_length = 151
    resources:
        mem_mb = 80000
    shell: """
        scripts/classify-kmers -jf {input.jf} -op {output.jf_position} -oq {output.jf_quality} -or {output.jf_read} \
        -reads {input.ss_fastq} -length-to-track-5p {params.read_length} -length-to-track-3p 0
    """

rule calculate_kmer_error_rate:
    input:
        jf_read = "results/kmer/{id}/read{read}/{id}_read{read}.tsv"
    output:
        kmer_accuarcy = temp("results/kmer/{id}/tmp/kmer_accuracy_{read}.txt")
    shell: """
        true=$(zcat {input.jf_read} | tail -n +2 | cut -f 16 | paste -sd+ | bc)
        total=$(zcat {input.jf_read} | tail -n +2 | cut -f 15 | paste -sd+ | bc)
        err_prop=$(echo "scale=4; (1-$true/$total)" | bc)
        formatted_result=$(printf "%06.4f" $err_prop)
        echo -e "{wildcards.id}\t$formatted_result" >> {output.kmer_accuarcy}
    """

rule aggregate_kmer_error_rate_1:
    input:
        files = expand("results/kmer/{id}/tmp/kmer_accuracy_1.txt", id = ids_1x_all)
    output:
        kmer_accuracy = "results/kmer/kmer_accuracy_read1.txt"
    shell: """
        cat {input.files} >> {output.kmer_accuarcy}
    """

rule aggregate_kmer_error_rate_2:
    input:
        files = expand("results/kmer/{id}/tmp/kmer_accuracy_2.txt", id = ids_1x_all)
    output:
        kmer_accuracy = "results/kmer/kmer_accuracy_read2.txt"
    shell: """
        cat {input.files} >> {output.kmer_accuarcy}
    """

rule plot_kmer_position:
    input:
        code = "scripts/plot_kmer_position.py",
        jf_pos1 = "results/kmer/{id}/read1/{id}_position1.tsv",
        jf_pos2 = "results/kmer/{id}/read2/{id}_position2.tsv"
    output:
        graph_kmer_position = "graphs/kmer_position/{id}_kmer_position.png"
    shell: """
	python {input.code} {wildcards.id}
    """

rule get_fragment_length:
    input:
        ss_bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        fragment_length = "results/kmer/{id}/fragment_length.tsv"
    shell: """
        samtools view -h {input.ss_bam} | grep -v '^@' | cut -f1,9 | \
        awk '{$2 = ($2 < 0 ? -$2 : $2); $2 = int(($2 + 49) / 50) * 50; print}' | \
        tr ' ' '\t' | \
        awk '$2 != 0' | \
        awk '!seen[$1]++' > {output.fragment_length}
    """

rule calculate_per_bin_kmer_error_rate:
    input:
        fragment_length = "results/kmer/{id}/fragment_length.tsv",
        reads = "results/kmer/{id}/read{read}/{id}_read{read}.tsv",
        script = "scripts/calculate_per_bin_kmer_error_rate.py"
    output:
        per_bin_kmer_accuarcy = "results/kmer/{id}/read{read}/per_bin_kmer_error_rate_read{read}.txt"
    resources: mem_mb = 50000
    shell: """
        python {input.script} {wildcards.id} {wildcards.read} {input.reads}
    """

rule aggregate_per_bin_kmer_error_rate_1:
    input:
        files = expand("results/kmer/{id}/read1/per_bin_kmer_error_rate_read1.txt", id = ids_1x_all)
    output:
        per_bin_kmer_accuracy1 = "results/kmer/per_bin_kmer_accuracy_read1.txt"
    shell: """
        cat {input.files} >> {output.per_bin_kmer_accuracy1}
    """

rule aggregate_per_bin_kmer_error_rate_2:
    input:
        files = expand("results/kmer/{id}/read2/per_bin_kmer_error_rate_read2.txt", id = ids_1x_all)
    output:
        per_bin_kmer_accuracy2 = "results/kmer/per_bin_kmer_accuracy_read2.txt"
    shell: """
        cat {input.files} >> {output.per_bin_kmer_accuracy2}
    """

