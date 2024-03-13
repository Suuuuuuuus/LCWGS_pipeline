configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
test_hc = samples_lc[:2]

hc_panel = config["hc_panel"]

rule extract_hla_alleles:
    input:
        vcf = "results/chip/imputed/20240313_no_qc/chr6.dose.vcf.gz"
    output:
        hla_all = temp("results/chip/imputed/20240313_no_qc/hla.vcf"),
        hla_all_gz = "results/chip/imputed/20240313_no_qc/hla.vcf.gz",
        two_field = temp("results/chip/imputed/20240313_no_qc/hla.two.vcf"),
        two_field_gz = "results/chip/imputed/20240313_no_qc/hla.two.vcf.gz"
    resources: mem = '10G'
    shell: """
        zgrep '^#' {input.vcf} > {output.hla_all}
        zgrep '^#' {input.vcf} > {output.two_field}

        zgrep -v '^#' {input.vcf} | awk '$3 ~ /HLA/' >> {output.hla_all} 
        zgrep -v '^#' {input.vcf} | awk '$3 ~ /HLA/ && split($3, a, ":") == 2' >> {output.two_field}

        gzip {output.hla_all}
        gzip {output.two_field}
    """

rule convert_to_manifest:
    input:
        bpm = "data/manifest/GDA-8v1-0_A1.bpm",
        egt = "data/manifest/GDA-8v1-0_A1_ClusterFile.egt"
    output:
        manifest = "data/manifest/manifest.csv"
    resources: mem = '10G'
    params:
        picard = tools['picard']
    shell: """
        {params.picard} BpmToNormalizationManifestCsv -I {input.bpm} -CF {input.egt} -O {output.manifest}
    """


 