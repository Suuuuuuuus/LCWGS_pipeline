configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"
# include: "hla_imputation_wip.smk"

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
        ref_panel = "data/tmp/a.tmp1.bam"

rule HLA_realignment:
    input:
        bam = "data/tmp/IDT0555.bam"
    output:
        bam = "data/tmp/a.bam",
        fastq1 = temp("data/tmp/a_1.fastq"),
        fastq2 = temp("data/tmp/a_2.fastq"),
        readlist = temp("data/tmp/a_readlist.txt"),
        tmp1 = "data/tmp/a.tmp1.bam"
    resources:
        mem = '40G'
    params: 
        picard = tools["picard"]
    threads: 8
    shell: """
        samtools view -h {input.bam} "chr6:25000000-34000000" | \
        samtools sort -n - -o {output.bam}
        
        samtools view {output.bam} | cut -f 1 | uniq > {output.readlist} 

        {params.picard} FilterSamReads \
        -I {input.bam} \
        -O {output.bam} \
        -READ_LIST_FILE {output.readlist} \
        -FILTER includeReadList

        samtools sort -n {output.bam} -o {output.tmp1}
        bedtools bamtofastq -i {output.tmp1} -fq {output.fastq1} -fq2 {output.fastq2}
    """

# rule HLA_realignment:
#     input:
#         bam = "data/bams/{id}.bam",
#         reference = "/well/band/users/rbx225/recyclable_files/hla_reference_files/fasta/v{IPD_IMGT_version}_{panel}/HLA.fasta"
#     output:
#         bam = "data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}.bam",
#         fastq1 = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}_1.fastq"),
#         fastq2 = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}_2.fastq"),
#         readlist = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}_readlist.txt"),
#         tmp1 = temp("data/realigned_bams/v{IPD_IMGT_version}_{panel}/{id}.tmp1.bam")
#     resources:
#         mem = '40G'
#     params: 
#         sample = "{id}",
#         picard = tools["picard"]
#     threads: 8
#     shell: """
#         mkdir -p data/realigned_bams/v{wildcards.IPD_IMGT_version}_{wildcards.panel}

#         samtools view -h {input.bam} "chr6:25000000-34000000" | \
#         samtools sort -n - -o {output.bam}
        
#         samtools view {output.bam} | cut -f 1 | uniq > {output.readlist} 

#         {params.picard} FilterSamReads \
#         -I {output.bam} \
#         -O {output.tmp1} \
#         -READ_LIST_FILE {output.readlist} \
#         -FILTER includeReadList

#         bedtools bamtofastq -i {output.tmp1} -fq {output.fastq1} -fq2 {output.fastq2}

#         bwa mem -t {threads} -a {input.reference} {output.fastq1} {output.fastq2} | \
#         samtools view -b -o {output.tmp1}
        
#         samtools sort -@{threads} -m 1G -o {output.bam} {output.tmp1}
#         samtools index {output.bam}
#     """