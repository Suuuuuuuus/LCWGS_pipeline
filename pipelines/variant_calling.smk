configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import os
import numpy as np
import sys
sys.path.append("scripts")
import lcwgSus

samples_hc = list(pd.read_table(config['samples_hc'], header = None, names = ['Code'])['Code'].values)
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
seq_names = list(sample_linker['Seq_Name'].values)
chip_names = list(sample_linker['Chip_Name'].values)
sample_names = list(sample_linker['Sample_Name'].values)

chromosomes = [i for i in range(1,23)]

test_hc = ids_1x_all[:2]

rule GATK_prepare_reference:
    input:
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"],
        bqsr_known_sites = config["bqsr_known_sites"]
    output:
        fai = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.fai" if concatenate else "data/references/GRCh38.fa.fai",
        dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict" if concatenate else "data/references/GRCh38.dict",
        bqsr_known_sites = [file + ".tbi" for file in config["bqsr_known_sites"]]
    resources: mem = '10G'
    shell: """
        samtools faidx {input.reference}
        picard CreateSequenceDictionary \
        R={input.reference} \
        O={output.dict}

        for i in {input.bqsr_known_sites}; do
            gatk IndexFeatureFile -I $i
        done
    """

rule get_bqsr_report:
    input:
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        bqsr_report = "results/call/BQSR/BQSR_reports/{hc}.chr{chr}.BQSR.report"
    params:
        bqsr_known_sites = config["bqsr_known_sites"]
    run:
        cmd = ""
        for file in params.bqsr_known_sites:
            cmd = cmd + "--known-sites " + file + " "
        shell("""
            gatk --java-options "-Xmx8G" BaseRecalibrator \
            -I {dedup_bam_chunk} \
            -R {ref} \
            -O {report} \
            {cmd}
        """.format(cmd = cmd, dedup_bam_chunk = input.dedup_bam_chunk, ref = input.reference, report = output.bqsr_report))

rule apply_bqsr:
    input:
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        bqsr_report = rules.get_bqsr_report.output.bqsr_report,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        recal_bam = "data/recal_bams/{hc}.chr{chr}.recal.bam",
        recal_bai = "data/recal_bams/{hc}.chr{chr}.recal.bam.bai"
    shell: """
        gatk --java-options "-Xmx8G" ApplyBQSR \
        -I {input.dedup_bam_chunk} \
        -R {input.reference} \
        --bqsr-recal-file {input.bqsr_report} \
        -O {output.recal_bam}

        samtools index {output.recal_bam}
    """

rule haplotype_call:
    input:
        recal_bam = rules.apply_bqsr.output.recal_bam,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        gvcf = "results/call/vcfs/regions/{hc}/{hc}.chr{chr}.gvcf.vcf.gz"
    resources: mem = '20G'
    threads: 8
    shell: """
        gatk --java-options "-Xmx20G" HaplotypeCaller \
        -R {input.reference} \
        -I {input.recal_bam} \
        -O {output.gvcf} \
        -ERC GVCF
    """

rule genomics_db_import:
    input:
        gvcfs = expand("results/call/vcfs/regions/{hc}/{hc}.chr{chr}.gvcf.vcf.gz", hc = test_hc, allow_missing = True),
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        temp(directory("results/call/tmp/{chr}.combined.db"))
    resources:
        mem = '15G'
    threads: 4
    run:
        gvcf_files = ""
        for gvcf in input.gvcfs:
            gvcf_files = gvcf_files + "--variant " + gvcf + " "
        shell("""
        gatk --java-options "-Xmx55g -Xms2g" GenomicsDBImport \
        -R {input_ref} \
        {gvcfs} \
        --tmp-dir=/well/band/users/rbx225/GAMCC/results/call/tmp/ \
        --genomicsdb-workspace-path {output}
        """.format(gvcfs = gvcf_files, input_ref = input.reference, output = output))

rule genotype_gvcf:
    input:
        gvcfs = rules.genomics_db_import.output,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        called = "results/call/vcfs/regions/chr{chr}.vcf.gz",
        index = "results/call/vcfs/regions/chr{chr}.vcf.gz.tbi"
    resources:
        mem = '15G'
    threads: 4
    shell: """
        gatk --java-options "-Xmx55g" GenotypeGVCFs  \
        -R {input.reference} \
        -V gendb://{input.gvcfs} \
        -O {output.called}

        bcftools index -t {output.called}
    """

'''
rule concat_hc_vcfs:
    input:
        vcfs = expand("results/call/vcfs/regions/chr{chr}.vcf.gz", chr = chromosomes),
        tbis = expand("results/call/vcfs/regions/chr{chr}.vcf.gz.tbi", chr = chromosomes)
    output:
        vcf = "results/call/vcfs/{hc}.vcf.gz" # Need to expand this if we're calling at different sites
    threads: 4
    resources: mem = '10G'
    shell: """
        bcftools concat --ligate -Oz -o {output.vcf} {input.vcfs}
    """
'''
