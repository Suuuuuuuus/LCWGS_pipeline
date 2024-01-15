configfile: "pipelines/config.json"

# Spliting fastq files
rule split_fastq:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        # dirs = directory("data/fastq/tmp/{id}/")，
        flag = temp("data/fastq/tmp/{id}/flag.txt")
    threads: 1
    params:
        chunk_size = config["fastq_chunk_size"]
    shell: """
        seqkit split2 -1 {input.fastq1} -2 {input.fastq2} -s {params.chunk_size} -O "data/fastq/tmp/{wildcards.id}" -f -e .gz
        echo "done!" > {output.flag}
    """

rule make_fastq_tsv:
    input:
        flag = "data/fastq/tmp/{id}/flag.txt"
    output:
        fastq_lsts = "data/file_lsts/hc_fastq_split/{id}_split.txt"
    threads: 1
    params:
        tmpdir = "data/fastq/tmp/{id}/"
    shell: """
        ls {params.tmpdir} | sed 's/{params.tmpdir}/g' | sed 's/.fastq.gz//g' > {output.fastq_lsts}
    """