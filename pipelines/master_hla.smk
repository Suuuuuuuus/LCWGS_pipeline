include: "hla.smk"
#include: "alignment.smk"
#include: "post_hla.smk"
include: "hla_ref_panel.smk"
include: "hla_imputation_wip.smk"
include: "phasing.smk"

include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append("/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/")
import lcwgsus
from lcwgsus.variables import *
from hla_phase import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
studies = ['1KG', 'GAMCC']

rule hla_imputation_prep_all:
    input:
        bamlist = expand("results/hla/imputation/bamlists/bamlist{num}.txt", num = bam_numbers),
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_v{IPD_IMGT_version}/HLA{gene}fullallelesfilledin.RData", gene = hla_genes, IPD_IMGT_version = IPD_IMGT_versions)

rule hla_imputation_all:
    input:
        hla_imputed = expand("results/hla/imputation/batches/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers)

two_stage_hla_vcf_outdir = config["two_stage_hla_vcf_outdir"]

rule post_hla_all:
    input:
        lifted = "results/hla/reference/multiEth_sites.b38.vcf.gz",
        two_stage_vcf = expand("{two_stage_hla_vcf_outdir}chr6.vcf.gz", two_stage_hla_vcf_outdir = two_stage_hla_vcf_outdir)

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

rule phasing_all:
    input:
        phased_vcf = "results/phasing/HLA_1KG_BEAGLE/phased.1KG.chr6.vcf.gz"

rule hla_ref_panel_all:
    input:
        # hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        # legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        # samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples",
        # RData = "results/hla_tests/quilt.hrc.hla.all.haplotypes.RData",
        # bam_all = "results/hla_tests/bamlist.txt",
        # html = expand("results/phasing/html/{study}-{gene}.html", gene = HLA_GENES, study = studies),
        phase_df = expand("results/phasing/phased_dfs/{study}-{gene}.tsv", gene = HLA_GENES, study = studies),

        fv_vcf = expand("results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf.gz", chr = chromosome),
        mGen_chunk_vcfs = [mGen_chunk_vcf_lst],
        oneKG_gamcc = expand("results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz", chr = chromosome)

rule hla_imputation_wip_all:
    input:
        bamlist = expand("results/hla/imputation/bamlists_fv/bamlist{num}.txt", num = bam_numbers),
        merged_ref_hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        merged_ref_legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        merged_ref_sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples",

        ref_panel_db = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
#        ref_panel_merged_ref = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, id = samples_fv),
        ref_panel_method = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
#        ref_panel_optimal = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, id = samples_fv),

        imputed_db = expand("results/hla/imputation/QUILT_HLA_result_db/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers),
#        imputed_merged_ref = expand("results/hla/imputation/QUILT_HLA_result_merged_ref/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        imputed_method = expand("results/hla/imputation/QUILT_HLA_result_method/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers),
#        imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),

