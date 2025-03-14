configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')

rule all:
    input:
        bam = expand("data/realigned_bams/v3390_oneKG/{id}.bam", id = ['IDT0555']),
        fastq1 = expand("data/realigned_bams/v3390_oneKG/{id}_1.fastq", id = ['IDT0555']),
        fastq2 = expand("data/realigned_bams/v3390_oneKG/{id}_2.fastq", id = ['IDT0555']),
        readlist = expand("data/realigned_bams/v3390_oneKG/{id}_readlist.txt", id = ['IDT0555']),
        tmp1 = expand("data/realigned_bams/v3390_oneKG/{id}.tmp1.bam", id = ['IDT0555'])

rule HLA_realignment:
    input:
        bam = "data/bams/{id}.bam",
        reference = "/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v3390_oneKG/HLA.fasta"
    output:
        bam = "data/realigned_bams/v3390_oneKG/{id}.bam",
        fastq1 = "data/realigned_bams/v3390_oneKG/{id}_1.fastq",
        fastq2 = "data/realigned_bams/v3390_oneKG/{id}_2.fastq",
        readlist = "data/realigned_bams/v3390_oneKG/{id}_readlist.txt",
        tmp1 = "data/realigned_bams/v3390_oneKG/{id}.tmp1.bam"
    resources:
        mem = '40G'
    params: 
        picard = tools["picard"]
    threads: 4
    shell: """
        mkdir -p data/realigned_bams/v3390_oneKG/

        samtools view {input.bam} "chr6:25000000-34000000" | \
        cut -f 1 | uniq > {output.readlist} 

        {params.picard} FilterSamReads \
        -I {input.bam} \
        -O {output.bam} \
        -READ_LIST_FILE {output.readlist} \
        -FILTER includeReadList

        samtools view -F 2308 -h {output.bam} | \
        samtools sort -n - -o {output.tmp1}

        bedtools bamtofastq -i {output.tmp1} -fq {output.fastq1} -fq2 {output.fastq2}

        bwa mem -t {threads} -T 0 -h 10000 10000 -a {input.reference} {output.fastq1} {output.fastq2} | \
        samtools view -b -o {output.tmp1}
        
        samtools sort -@{threads} -m 1G -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """