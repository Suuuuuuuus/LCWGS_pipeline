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
# samples_lc = read_tsv_as_lst(config['samples_lc'])
# test_hc = samples_lc[:2]
chromosome = [i for i in range(1,23)]
variant_types = ['snps', 'indels']

rule GATK_prepare_reference:
    input:
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"],
        bqsr_known_sites = config["bqsr_known_sites"],
        gatk_to_index = config["gatk_to_index"]
    output:
        fai = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.fai" if concatenate else "data/references/GRCh38.fa.fai",
        dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict" if concatenate else "data/references/GRCh38.dict",
        bqsr_known_sites = [file + ".tbi" for file in config["bqsr_known_sites"]],
        gatk_to_index = [file + ".tbi" for file in config["gatk_to_index"]]
    resources: mem = '10G'
    shell: """
        samtools faidx {input.reference}
        picard CreateSequenceDictionary \
        R={input.reference} \
        O={output.dict}

        for i in {input.bqsr_known_sites}; do
            gatk IndexFeatureFile -I $i
        done

        for i in {input.gatk_to_index}; do
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
        -OBI false \
        -O {output.recal_bam}

        samtools index {output.recal_bam}
    """

rule prepare_hc_bamlist:
    input:
       # bams = expand("data/recal_bams/{hc}.recal.bam", hc = samples_hc)
    output:
        bamlist = "results/call/bam.list"
    shell: """
        mkdir -p results/call
        ls data/recal_bams/*.recal.bam > {output.bamlist}
    """

rule GATK_chunk_reference:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        ref_vcf = f"data/ref_panel/{hc_panel}/{hc_panel}.chr{{chr}}.vcf.gz" # Move this to config so that we can test sites at different ref panels
    output:
        empty_vcf1 = temp("results/call/tmp/ref/empty1_{type}_{regionStart}.{regionEnd}_chr{chr}.vcf.gz"),
        empty_vcf2 = temp("results/call/tmp/ref/empty2_{type}_{regionStart}.{regionEnd}_chr{chr}.vcf.gz")
    resources: mem = '10G'
    shell: """
        mkdir -p results/call/tmp/ref/
        mkdir -p results/call/vcfs/{hc_panel}/
        file=$(head -n 1 {input.bamlist})

        bcftools view -G -v {wildcards.type} -r chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} -Oz -o {output.empty_vcf1} {input.ref_vcf}
        
        gatk IndexFeatureFile -I {output.empty_vcf1}

        gatk UpdateVCFSequenceDictionary \
        -V {output.empty_vcf1} \
        --source-dictionary $file \
        --output {output.empty_vcf2} \
        --replace true
    """

rule haplotype_call:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        empty_vcf2 = rules.GATK_chunk_reference.output.empty_vcf2
    output:
        vcf = f"results/call/vcfs/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources: mem = '20G'
    params:
        padding = 300
    threads: 4
    shell: """
        if [[ {wildcards.type} == "indels" ]]
        then
            gatk --java-options "-Xmx20G -Xms20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {input.empty_vcf2} \
            --alleles {input.empty_vcf2} \
            --native-pair-hmm-threads 4 \
            --assembly-region-padding {params.padding} \
            --output-mode EMIT_VARIANTS_ONLY
        else
            gatk --java-options "-Xmx20G -Xms20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {input.empty_vcf2} \
            --alleles {input.empty_vcf2} \
            --native-pair-hmm-threads 4 \
            --output-mode EMIT_VARIANTS_ONLY
        fi

        rm "{input.empty_vcf2}.tbi"
    """

rule get_vqsr_report:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.haplotype_call.output.vcf
    output:
        tranch = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.{{regionStart}}.{{regionEnd}}.tranch",
        recal = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.{{regionStart}}.{{regionEnd}}.recal"
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
            -mode SNP -O {output.recal} --tranches-file {output.tranch}
        else
            gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.oneKG_indels} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL -O {output.recal} --tranches-file {output.tranch}
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
        recal_vcf = f"results/call/recal_vcf/{hc_panel}/regions/{hc_panel}.{{type}}.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/call/recal_vcf/{hc_panel}/regions/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.recal_vcf}
        else
            gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.recal_vcf}
        fi
    """

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if os.path.exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file)

vcfs_to_concat={}
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    vcfs_to_concat[str(chr)] = {}
    for t in variant_types:
        file_ary = []
        for i in range(0, start.__len__()):
            regionStart=start[i]
            regionEnd=end[i]
            file="results/call/recal_vcf/" + hc_panel + "/regions/" + hc_panel + "." + t +  ".chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            file_ary.append(file)
        vcfs_to_concat[str(chr)][t] = file_ary

rule concat_hc_vcfs:
    input:
        vcfs = lambda wildcards: vcfs_to_concat[str(wildcards.chr)][str(wildcards.type)]
    output:
        vcf = f"results/call/recal_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    threads: 4
    resources: mem = '10G'
    shell: """
        bcftools concat --ligate -Oz -o {output.vcf} {input.vcfs}
    """
