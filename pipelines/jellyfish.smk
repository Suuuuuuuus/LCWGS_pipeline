configfile: "pipelines/config.json"

import pandas as pd
config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)

rule classify_kmers:
    input:
        fastq = "data/subsampled_fastq/{id}_subsampled_{read}.fastq",
        jf = "data/jellyfish/NA12878_31mer.jf"
    output:
        jf_read = "results/kmer/{id}/read{read}/{id}_read{read}.tsv.gz",
        jf_quality = "results/kmer/{id}/read{read}/{id}_quality{read}.tsv.gz",
        jf_position = "results/kmer/{id}/read{read}/{id}_position{read}.tsv.gz"
    params:
        read_length = 151
    resources:
        mem_mb = 80000
    shell: """
        scripts/classify-kmers -jf {input.jf} -op {output.jf_position} -oq {output.jf_quality} -or {output.jf_read} \
        -reads {input.fastq} -length-to-track-5p {params.read_length} -length-to-track-3p 0
    """

rule calculate_kmer_error_rate:
    input:
        jf_read = "results/kmer/{id}/read{read}/{id}_read{read}.tsv.gz"
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
        jf_pos1 = "results/kmer/{id}/read1/{id}_position1.tsv.gz",
        jf_pos2 = "results/kmer/{id}/read2/{id}_position2.tsv.gz"
    output:
        graph_kmer_position = "graphs/kmer_position/{id}_kmer_position.png"
    shell: """
	python {input.code} {wildcards.id}
    """
