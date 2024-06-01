include: "hla.smk"
include: "alignment.smk"
include: "post_hla.smk"

include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']

hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_ref_files/"

rule hla_imputation_prep_all:
    input:
        bamlist = "results/hla/imputation/bamlist.txt",
        # ref_panel = expand(f"{hla_ref_panel_outdir}HLA{{gene}}fullallelesfilledin.RData", gene = hla_genes)

# rule hla_imputation_prep_alt_all:
#     input:
#         chrs = expand("data/hla_bams_alt/{id}.chr6.bam", id = samples_lc),
#         # bamlist = "results/hla/imputation/bamlist.txt",
#         aligned_bams = expand("data/hla_bams_alt/{id}.chr6.tmp.bam", id = samples_lc),
#         # hla_imputed = expand("results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes)

rule hla_imputation_all:
    input:
        hla_imputed = expand("results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes)

# rule hla_imputation_alt_all:
#     input:
#         hla_imputed = expand("results/hla/imputation/genes/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes),
#         # bamlist = "results/hla/imputation/bamlist_alt.txt"

two_stage_vcf_outdir = config["two_stage_vcf_outdir"]
three_stage_vcf_outdir = config["three_stage_vcf_outdir"]

rule post_hla_all:
    input:
        two_stage_vcf = expand("{two_stage_vcf_outdir}chr6.vcf.gz", two_stage_vcf_outdir = two_stage_vcf_outdir),
        three_stage_vcf = expand("{three_stage_vcf_outdir}chr6.vcf.gz", three_stage_vcf_outdir = three_stage_vcf_outdir)