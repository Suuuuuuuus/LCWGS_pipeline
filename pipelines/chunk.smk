configfile: "pipelines/config.json"
include: "auxiliary.smk"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

chunks = read_tsv_as_lst("data/bedgraph/bam_chunks.bed")
    
# Spliting fastq files
rule split_fastq:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        dirs = directory("data/fastq/tmp/{id}/"),
        flag = temp("data/fastq/tmp/{id}/flag.txt")
    threads: 1
    params:
        chunk_size = config["fastq_chunk_size"]
    shell: """
        mkdir -p "data/fastq/tmp/{wildcards.id}/"
        seqkit split2 -1 {input.fastq1} -2 {input.fastq2} -s {params.chunk_size} -O "data/fastq/tmp/{wildcards.id}" -f -e .gz
        echo "done!" > {output.flag}
    """

rule make_fastq_tsv:
    input:
        flag = "data/fastq/tmp/{id}/flag.txt"
    output:
        fastq_lsts = "data/file_lsts/hc_fastq_split/{id}_split.tsv"
    threads: 1
    params:
        tmpdir = "data/fastq/tmp/{id}/"
    shell: """
        ls {params.tmpdir}*_1.part* | sed 's/_1.part//g' | sed 's/.fastq.gz/_1.fastq.gz/g' > "{params.tmpdir}tmp1.txt"
        ls {params.tmpdir}*_2.part* | sed 's/_2.part//g' | sed 's/.fastq.gz/_2.fastq.gz/g' > "{params.tmpdir}tmp2.txt"
        n=1
        for i in $(ls {params.tmpdir}*_1.part*); do
            mv $i $(sed -n "${{n}}p" "{params.tmpdir}tmp1.txt")
            n=$((n+1))
        done
        n=1
        for i in $(ls {params.tmpdir}*_2.part*);
            mv $i $(sed -n "${{n}}p" "{params.tmpdir}tmp2.txt")
            n=$((n+1))
        done
        cat "{params.tmpdir}tmp1.txt" | rev | cut -d '/' -f1 | rev | sed 's/_1.fastq.gz//g' > {output.fastq_lsts}
        rm "{params.tmpdir}tmp1.txt" "{params.tmpdir}tmp2.txt"
    """

rule split_bams:
    input:
        bam = "data/merge_bams/{id}.bam"
    output:
        bam_chunk = temp("data/chunk_bams/tmp/{id}/{id}.chr{chr}.bam")
    threads: 8
    resources: mem = '100G'
    shell: """
        mkdir -p data/chunk_bams/tmp/
        samtools view -h {input.bam} {wildcards.chr} | \
        samtools sort -n | \
        samtools fixmate -m | \
        samtools sort | \
        samtools markdup -O BAM -o {output.bam_chunk}
    """