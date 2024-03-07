configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import os
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

hc_panel = config["hc_panel"]
concatenate = config["concatenate"]

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
test_hc = samples_lc[:2]
chromosome = [i for i in range(1,23)]

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
        dedup_bam = "data/merge_bams/tmp/{hc}.bam",
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        bqsr_report = "results/call/BQSR/BQSR_reports/{hc}.BQSR.report"
    params:
        bqsr_known_sites = config["bqsr_known_sites"]
    resources:
        mem = '10G'
    run:
        cmd = ""
        for file in params.bqsr_known_sites:
            cmd = cmd + "--known-sites " + file + " "
        shell("""
            gatk --java-options "-Xmx8G" BaseRecalibrator \
            -I {dedup_bam} \
            -R {ref} \
            -O {report} \
            {cmd}
        """.format(cmd = cmd, dedup_bam = input.dedup_bam, ref = input.reference, report = output.bqsr_report))

rule apply_bqsr:
    input:
        dedup_bam = "data/merge_bams/tmp/{hc}.bam",
        bqsr_report = rules.get_bqsr_report.output.bqsr_report,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        recal_bam = "data/recal_bams/{hc}.recal.bam",
        recal_bai = "data/recal_bams/{hc}.recal.bam.bai"
    resources:
        mem = '10G'
    shell: """
        gatk --java-options "-Xmx8G" ApplyBQSR \
        -I {input.dedup_bam} \
        -R {input.reference} \
        --bqsr-recal-file {input.bqsr_report} \
        -O {output.recal_bam}

        samtools index {output.recal_bam}
    """

rule prepare_hc_bamlist:
    input:
        bams = expand("data/recal_bams/{hc}.recal.bam", hc = test_hc)
    output:
        bamlist = "results/call/bam.list"
    shell: """
        mkdir -p results/call
        ls data/recal_bams/*.recal.bam > {output.bamlist}
    """

rule haplotype_call:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        ref_vcf = f"data/ref_panel/{hc_panel}/{hc_panel}.chr{{chr}}.vcf.gz" # Move this to config so that we can test sites at different ref panels
    output:
        vcf = f"results/call/vcfs/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz",
        empty_vcf1 = temp("results/call/tmp/ref/empty1_{type}_chr{chr}.vcf.gz"),
        empty_vcf2 = temp("results/call/tmp/ref/empty2_{type}_chr{chr}.vcf.gz")
    resources: mem = '20G'
    threads: 8
    shell: """
        mkdir -p results/call/tmp/ref/
        mkdir -p results/call/vcfs/{hc_panel}/
        file=$(head -n 1 {input.bamlist})

        bcftools view -G {input.ref_vcf} | bcftools view -v {wildcards.type} -Oz -o {output.empty_vcf1}
        gatk IndexFeatureFile -I {output.empty_vcf1}

        gatk UpdateVCFSequenceDictionary \
        -V {output.empty_vcf1} \
        --source-dictionary $file \
        --output {output.empty_vcf2} \
        --replace true

        gatk --java-options "-Xmx20G" HaplotypeCaller \
        -R {input.reference} \
        -I {input.bamlist} \
        -O {output.vcf} \
        -L {output.empty_vcf2} \
        --alleles {output.empty_vcf2} \
        --output-mode EMIT_VARIANTS_ONLY

        rm "{output.empty_vcf1}.tbi" "{output.empty_vcf2}.tbi"
    """

rule get_vqsr_report:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.haplotype_call.output.vcf
    output:
        tranch = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.tranches",
        recal = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.recal"
    params:
        hapmap = "data/GATK_resource_bundle/hapmap_3.3.hg38.vcf.gz",
        omni = "data/GATK_resource_bundle/1000G_omni2.5.hg38.vcf.gz",
        oneKG_snps = "data/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        GRCh38_indels = "data/GATK_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz",
        oneKG_indels = "data/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "data/GATK_resource_bundle/dbsnp_146.hg38.vcf.gz"
    resources:
        mem = '30G'
    shell: """
        mkdir -p results/call/VQSR/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.oneKG_snps} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=7 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP -O {output.recal} --tranches-file {output.tranches}
        else
            gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.oneKG_indels} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL -O {output.recal} --tranches-file {output.tranches}
        fi
    """

rule apply_vqsr:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.haplotype_call.output.vcf,
        tranch = rules.get_vqsr_report.output.tranch,
        recal = rules.get_vqsr_report.output.recal
    output:
        recal_vcf = f"results/call/recal_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    params:
        mode = "SNP" if "{type}" == "snps" else "INDEL"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/call/recal_vcf/{hc_panel}/

        gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
        -V {input.vcf} \
        --recal-file {input.recal} \
        --tranches-file {input.tranch} \
        --truth-sensitivity-filter-level 99.0 \
        --create-output-variant-index true \
        -mode {params.mode} \
        -O {output.recal_vcf}
    """