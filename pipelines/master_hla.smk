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

bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]

rule hla_imputation_prep_all:
    input:
        bamlist = expand("results/hla/imputation/bamlists/bamlist{num}.txt", num = bam_numbers),
        # ref_panel = expand(f"{hla_ref_panel_outdir}HLA{{gene}}fullallelesfilledin.RData", gene = hla_genes)

rule hla_imputation_all:
    input:
        hla_imputed = expand("results/hla/imputation/batches/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers)

two_stage_hla_vcf_outdir = config["two_stage_hla_vcf_outdir"]

rule post_hla_all:
    input:
        lifted = "results/hla/reference/multiEth_sites.b38.vcf.gz",
        two_stage_vcf = expand("{two_stage_hla_vcf_outdir}chr6.vcf.gz", two_stage_hla_vcf_outdir = two_stage_hla_vcf_outdir)