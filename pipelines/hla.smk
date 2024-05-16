configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
QUILT_HOME = config["QUILT_HOME"]

rule index:
    input:
        reference = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        indexed = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb"
    resources:
        mem = '50G'
    threads: 8
    shell: """
        bwa index {input.reference}
    """

rule alignment:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz",
        reference = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        index = rules.index.output.indexed
    output:
        bam = temp("data/bams/tmp/{id}.bam")
    resources:
        mem = '10G'
    threads: 8
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
    """

rule hla_imputation_preprocess:
    input:
        bam = "data/bams/tmp/{id}.bam"
    output:
        tmp = temp("data/hla_bams/{id}.tmp.bam"),
        chr = temp("data/hla_bams/{id}.chr6.tmp.bam")
    params:
        verbosity = "ERROR",
        sample = "{id}"
    threads: 4
    resources: mem='30G'
    shell: """
        mkdir -p data/hla_bams/

        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {input.bam} \
        -O {output.tmp} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}

        samtools sort -@4 -m 1G -o {output.chr} {output.tmp}

        samtools index {output.chr}
    """

rule hla_clean_bam:
    input:
        bam = rules.hla_imputation_preprocess.output.chr
    output:
        bam = "data/hla_bams/{id}.chr6.bam",
        bai = "data/hla_bams/{id}.chr6.bam.bai",
        sam = temp("data/hla_bams/{id}.chr6.sam"),
        tmp1 = temp("data/hla/bams/{id}.tmp1.bam"),
        metric = temp("data/hla/bams/{id}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/hla/bams/tmp/{id}/",
        sample = "{id}"
    shell: """
        mkdir -p {params.tmpdir}

        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {input.bam}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

# Filter chr6 and HLA contigs as well as recode QUAL strs

        samtools view -H {output.bam} > {output.sam}
        samtools view {output.bam} | \
        awk -F '\t' '($3~/chr6/||$3~/HLA/){{print}}' | \
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1 ~ /^@/) {{
                print $0
            }} else {{
                new_qual = ""
                qual = $11
                for (i = 1; i <= length(qual); i++) {{
                    q = substr(qual, i, 1)
                    if (q ~ /[@ABCDEF]/) {{
                        q = "?"
                    }}
                    new_qual = new_qual q
                }}
                $11 = new_qual
                print $0
            }}
        }}' >> {output.sam}
        samtools view -bS {output.sam} > {output.bam}

        samtools index {output.bam}
    """

rule prepare_hla_bamlist:
    input:
        bams = expand("data/hla_bams/{id}.chr6.bam", id = samples_lc)
    output:
        bamlist = "results/hla/imputation/bamlist.txt"
    shell: """
        mkdir -p results/hla/imputation/

        ls data/hla_bams/*.bam > {output.bamlist}
    """

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_ref_files/"
hla_ref_panel_outdir_original = "results/hla/imputation/ref_panel/QUILT_ref_files_original/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']

bamlist_file = "results/hla/imputation/bamlist.txt"
bamlist_test_file = "results/hla/imputation/test.txt"

rule prepare_hla_reference_panel:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD_IGMT.zip",
        fasta = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        genetic_map = f"{hla_ref_panel_indir}YRI/YRI-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand(f"{hla_ref_panel_outdir}HLA{{gene}}fullallelesfilledin.RData", gene = hla_genes)
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = "/well/band/users/rbx225/software/QUILT/hla_ancillary_files/hlagenes.txt"
    shell: """
        {params.quilt_hla_prep} \
        --outputdir={hla_ref_panel_outdir} \
        --nGen=100 \
        --hla_types_panel={input.hla_types_panel} \
        --ipd_igmt_alignments_zip_file={input.ipd_igmt} \
        --ref_fasta={input.fasta} \
        --refseq_table_file={params.refseq} \
        --full_regionStart=25587319 \
        --full_regionEnd=33629686 \
        --buffer=500000 \
        --region_exclude_file={params.region_exclude_file} \
        --genetic_map_file={input.genetic_map} \
        --reference_haplotype_file={input.hap} \
        --reference_legend_file={input.legend} \
        --reference_sample_file={input.sample} \
        --hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
        --nCores=6
    """

rule hla_imputation:
    input:
        bamlist = bamlist_file,
        ref_dir = hla_ref_panel_outdir_original
    output:
        vcf = "results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '30G'
    threads: 4
    params:
        quilt_hla = tools['quilt_hla']
    shell: """
        mkdir -p results/hla/imputation/genes/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/genes/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={QUILT_HOME}hla_ancillary_files/GRCh38_full_analysis_set_plus_decoy_hla.dict
    """
