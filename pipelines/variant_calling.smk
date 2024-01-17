configfile: "pipelines/config.json"

from os.path import exists
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

rule get_bqsr_report:
    input:
        bam = "data/dedup_bams/{id}.bam" if config["rmdup"] else "data/bams/{id}.bam",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
    output:
        bqsr_report = "results/call/BQSR/BQSR_reports/{id}_BQSR.report"
    params:
        bqsr_known_sites = config["bqsr_known_sites"]
    shell: """
        cmd="gatk --java-options "-Xmx8G" BaseRecalibrator \
        -I {input.bam} \
        -R {input.reference} \
        -O {output.bqsr_report}"

        for file in "${params.bqsr_known_sites}[@]"; do
            cmd+=" --known-sites $file"
        done
        
        eval "$cmd"
    """

rule apply_bqsr:
    input:
        bam = "data/dedup_bams/{id}.bam" if config["rmdup"] else "data/bams/{id}.bam",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"],
        bqsr_report = "results/call/BQSR/BQSR_reports/{id}_BQSR.report"
    output:
        recal_bam = "data/recal_bams/{id}_recal.bam"
    shell: """
        gatk --java-options "-Xmx8G" ApplyBQSR \
        -I {input.bam} \
        -R {input.reference} \
        --bqsr-recal-file {input.bqsr_report} \
        -O {output.recal_bam}
    """

rule index_recal:
	input:
		recal_bam = rules.apply_bqsr.output.recal_bam
	output:
		bai = "data/recal_bams/{id}_recal.bam.bai"
	resources: mem = '50G'
	shell: """
		samtools index {input.recal_bam}
	"""

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file) ## python is dumb

vcfs_all_samples={}
vcfs_to_concat={}
final_vcfs=[]
for sample in samples_hc:
    for chr in chromosomes:
        start=REGIONS[str(chr)]["start"]
        end=REGIONS[str(chr)]["end"]
        per_chr_vcfs=[]
        for i in range(0, start.__len__()):
            regionStart=start[i]
            regionEnd=end[i]
            file="results/call/vcfs/" + sample + "regions/" + sample + ".chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            per_chr_vcfs.append(file)
        vcfs_to_concat[str(chr)]=per_chr_vcfs
        final_vcfs.append("results/call/vcfs/" + sample + "/" + sample + ".chr" + str(chr) + ".vcf.gz")
    vcfs_all_samples[sample] = final_vcfs