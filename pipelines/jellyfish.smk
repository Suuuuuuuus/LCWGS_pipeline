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
        mem_mb = 60000
    shell: """
        scripts/classify_kmers -jf {input.jf} -op {output.jf_position} -oq {output.jf_quality} -or {output.jf_read} \
        -reads {input.fastq} -length-to-track-5p {params.read_length} -length-to-track-3p 0
    """

rule calculate_kmer_error_rate:
    input:
        jf_read = "results/kmer/{id}/read{read}/{id}_read{read}.tsv.gz"
    output:
        kmer_accuarcy = temp("results/kmer/{id}/tmp/kmer_accuracy_{read}.txt")
    shell: """
        true=$(zcat {input.jf_read} | tail -n +2 | cut -d ' ' -f 16 | paste -sd+ | bc)
        total=$(zcat {input.jf_read} | tail -n +2 | cut -d ' ' -f 15 | paste -sd+ | bc)
        err_prop=$(echo "scale=4; (1-$true/$total)" | bc)
        formatted_result=$(printf "%06.4f" $err_prop)
        echo -e "$formatted_result" >> {output.kmer_accuarcy}
    """

rule calculate_avg_kmer_error_rate:
    input:
        res1 = "results/kmer/{id}/tmp/kmer_accuracy_1.txt",
        res2 = "results/kmer/{id}/tmp/kmer_accuracy_2.txt"
    output:
        kmer_accuarcy = temp("results/kmer/{id}/tmp/kmer_accuracy.txt")
    params:
        tmpdir = "results/dup_rate/tmp/"
    shell: """
        sum_err=$(cut -f2 -d ' ' results/kmer/{wildcards.id}/tmp/* | paste -sd+ | bc)
        avg_err=$(echo "scale=4; $sum_err/2" | bc)
        formatted_result=$(printf "%06.4f" $avg_err)
        echo -e "{wildcards.id}\t$formatted_result" >> {output.kmer_accuarcy}
    """

rule aggregate_kmer_error_rate:
    input:
        files = expand("results/kmer/{id}/tmp/kmer_accuracy.txt", id = ids_1x_all)
    output:
        kmer_accuarcy = "results/kmer/kmer_accuracy.txt"
    shell: """
        cat {input.files} >> {output.kmer_accuarcy}
    """

'''
rule plot_satools_dup_rate:
    input:
        script = "scripts/plot_samtools_duplication_rate.py"
    output:
        graph_samtools_dup_rate = "graphs/samtools_duplication_rate.png"
    shell: """
        python {input.script}
    """
'''