configfile: "pipelines/config.json"
#include: "chunk.smk"
#include: "preprocess.smk"
#include: "reference.smk"
#include: "alignment.smk"

#include: "merge.smk"

include: "variant_calling.smk"

#include: "test.smk"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
test_hc = samples_lc[:2]

hc_panel = config["hc_panel"]

chromosome = [i for i in range(1,23)]

# The followings are global parameters:
clean_fastq = config['clean_fastq']
reheader = config['reheader']
concatenate = config['concatenate']

# chunks = read_tsv_as_lst("data/bedgraph/bam_chunks.bed")

rule chunk_all:
    input:
        # flag = expand("data/fastq/tmp/{hc}/flag.txt", hc = samples_hc),
        fastq_lsts = expand("data/file_lsts/hc_fastq_split/{hc}_split.tsv", hc = samples_hc)

samples_hc_split = []
for i in samples_hc:
    samples_hc_split = samples_hc_split + read_tsv_as_lst("data/file_lsts/hc_fastq_split/" + i + "_split.tsv")

rule preprocess_all:
    input:
        fwd_pair = expand("data/fastq_cleaned/{id}_1.fastq.gz", id = samples_hc_split),
        rev_pair = expand("data/fastq_cleaned/{id}_2.fastq.gz", id = samples_hc_split),
        fwd_unpair = expand("data/fastq_cleaned/{id}_unpaired_1.fastq.gz", id = samples_hc_split),
        rev_unpair = expand("data/fastq_cleaned/{id}_unpaired_2.fastq.gz", id = samples_hc_split)

rule reference_all:
    input:
        amb = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.amb" if concatenate else "data/references/GRCh38.fa.amb",
        ann = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.ann" if concatenate else "data/references/GRCh38.fa.ann",
        bwt = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.bwt" if concatenate else "data/references/GRCh38.fa.bwt",
        pac = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.pac" if concatenate else "data/references/GRCh38.fa.pac",
        sa = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.sa" if concatenate else "data/references/GRCh38.fa.sa"

rule fastqc_all:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html", id = samples_hc),
        html2 = expand("results/fastqc/{id}_2_fastqc.html", id = samples_hc),
        zip1 = expand("results/fastqc/{id}_1_fastqc.zip", id = samples_hc),
        zip2 = expand("results/fastqc/{id}_2_fastqc.zip", id = samples_hc),
        multiqc_hc = "results/fastqc/multiqc_hc/multiqc_report.html"

rule alignment_all:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_hc_split),
        bais = expand("data/bams/{id}.bam.bai", id = samples_hc_split)

rule merge_all:
    input:
        bams = expand("data/merge_bams/tmp/{hc}.bam", hc = samples_hc[:4]),
        bais = expand("data/merge_bams/tmp/{hc}.bam.bai", hc = samples_hc[:4])

variant_types = ['snps', 'indels']

rule variant_calling_all:
    input:
        fai = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.fai" if concatenate else "data/references/GRCh38.fa.fai",
        dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict" if concatenate else "data/references/GRCh38.dict",
        bqsr_known_sites = [file + ".tbi" for file in config["bqsr_known_sites"]],
        bqsr_reports = expand("results/call/BQSR/BQSR_reports/{hc}.BQSR.report", hc = test_hc),
        recal_bams = expand("data/recal_bams/{hc}.recal.bam", hc = test_hc),
        recal_bais = expand("data/recal_bams/{hc}.recal.bam.bai", hc = test_hc),
        bamlist = "results/call/bam.list",
        vcf = expand(f"results/call/vcfs/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz", chr = chromosome, type = variant_types),
        tranch = expand(f"results/call/VQSR/{hc_panel}/{{type}}.tranches", type = variant_types),
        recal = expand(f"results/call/VQSR/{hc_panel}/{{type}}.recal", type = variant_types),
        recal_vcf = expand(f"results/call/recal_vcf/{hc_panel}/{{type}}.vcf.gz", type = variant_types)
        
rule test_all:
    input:
        vcf = "results/tmp/{id}.{chr}.txt"
        # bam = expand("data/bams/tmp/{hc}.bam", hc = test_hc),
        # bai = expand("data/bams/tmp/{hc}.bam.bai", hc = test_hc)