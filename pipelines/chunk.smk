configfile: "pipelines/config.json"
include: "reference.smk"
include: "auxiliary.smk"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

# chunks = read_tsv_as_lst("data/bedgraph/bam_chunks.bed")
chromosome = [i for i in range(1,23)]
    
# Spliting fastq files
rule split_fastq:
    input:
        fastq1 = "data/fastq/{hc}_1.fastq.gz",
        fastq2 = "data/fastq/{hc}_2.fastq.gz"
    output:
        flag = touch("data/fastq/tmp/{hc}/flag.txt")
    threads: 1
    params:
        chunk_size = config["fastq_chunk_size"]
    shell: """
        mkdir -p "data/fastq/tmp/{wildcards.hc}/"
        seqkit split2 -1 {input.fastq1} -2 {input.fastq2} -s {params.chunk_size} -O "data/fastq/tmp/{wildcards.hc}" -f -e .gz
    """

rule make_fastq_tsv:
    input:
        flag = rules.split_fastq.output.flag
    output:
        fastq_lsts = "data/file_lsts/hc_fastq_split/{hc}_split.tsv"
    threads: 1
    params:
        tmpdir = "data/fastq/tmp/{hc}/"
    shell: """
        ls {params.tmpdir}*_1.part* | sed 's/_1.part//g' | sed 's/.fastq.gz/_1.fastq.gz/g' > "{params.tmpdir}tmp1.txt"
        ls {params.tmpdir}*_2.part* | sed 's/_2.part//g' | sed 's/.fastq.gz/_2.fastq.gz/g' > "{params.tmpdir}tmp2.txt"
        n=1
        for i in $(ls {params.tmpdir}*_1.part*); do
            mv $i $(sed -n "${{n}}p" "{params.tmpdir}tmp1.txt")
            n=$((n+1))
        done
        n=1
        for i in $(ls {params.tmpdir}*_2.part*); do
            mv $i $(sed -n "${{n}}p" "{params.tmpdir}tmp2.txt")
            n=$((n+1))
        done
        cat "{params.tmpdir}tmp1.txt" | rev | cut -d '/' -f1 | rev | sed 's/_1.fastq.gz//g' > {output.fastq_lsts}
        rm "{params.tmpdir}tmp1.txt" "{params.tmpdir}tmp2.txt"
    """