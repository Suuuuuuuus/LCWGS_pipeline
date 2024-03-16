configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
hc_panel = config["hc_panel"]
chromosome = [i for i in range(1,23)]
variant_types = ['snps', 'indels']
concatenate = config["concatenate"]

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
            file = "results/call/vcfs/" + hc_panel + "/" + hc_panel + "." + t +  ".chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            file_ary.append(file)
        vcfs_to_concat[str(chr)][t] = file_ary

rule concat_hc_vcfs:
    input:
        vcfs = lambda wildcards: vcfs_to_concat[str(wildcards.chr)][str(wildcards.type)]
    output:
        vcf = f"results/call/merge_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    threads: 4
    resources: mem = '20G'
    shell: """
        mkdir -p results/call/merge_vcf/{hc_panel}/

        bcftools concat -Oz -o {output.vcf} -a -d {wildcards.type} {input.vcfs}
        tabix {output.vcf}
    """

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

rule get_vqsr_report:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.concat_hc_vcfs.output.vcf
    output:
        tranch = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.tranch",
        recal = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.recal"
    params:
        hapmap = "data/GATK_resource_bundle/hapmap_3.3.hg38.vcf.gz",
        omni = "data/GATK_resource_bundle/1000G_omni2.5.hg38.vcf.gz",
        oneKG_snps = "data/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        GRCh38_indels = "data/GATK_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz",
        oneKG_indels = "data/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "data/GATK_resource_bundle/dbsnp_146.hg38.vcf.gz"
    resources:
        mem = '20G'
    shell: """
        mkdir -p results/call/VQSR/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \
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
            gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \
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
        vcf = rules.concat_hc_vcfs.output.vcf,
        tranch = rules.get_vqsr_report.output.tranch,
        recal = rules.get_vqsr_report.output.recal
    output:
        recal_vcf = f"results/call/recal_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    resources:
        mem = '20G'
    shell: """
        mkdir -p results/call/recal_vcf/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xmx20g -Xms20g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.recal_vcf}
        else
            gatk --java-options "-Xmx20g -Xms20g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.recal_vcf}
        fi
    """