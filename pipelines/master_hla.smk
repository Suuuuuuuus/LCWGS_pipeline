include: "hla.smk"
#include: "alignment.smk"
#include: "post_hla.smk"
include: "hla_ref_panel.smk"
#include: "hla_imputation_wip.smk"
include: "hla_imputation_method.smk"
#include: "hla_imputation_prep.smk"
include: "phasing.smk"
include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import io
import os
import re
import json
import pandas as pd
import numpy as np
import math
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/')
import lcwgsus

from lcwgsus.variables import *
from hla_phase import *
from hla_align_functions import *
from hla_align import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gm = read_tsv_as_lst('data/sample_tsvs/fv_gm_names.tsv')
samples_oneKG = read_tsv_as_lst("/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/samples_to_phase.tsv")
chromosome = [i for i in range(1,23)]
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']
# IPD_IMGT_versions = ['3390']
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
studies = ['1KG', 'GAMCC']
filters = ['strict', 'loose']
vcf_versions = ['30x', 'phase3_b38']
panels = ['oneKG', 'merged']

rule hla_all:
    input:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_v{IPD_IMGT_version}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, IPD_IMGT_version = IPD_IMGT_versions),
        imputed = expand("results/hla/imputation/QUILT_HLA_result_v{IPD_IMGT_version}/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", IPD_IMGT_version = IPD_IMGT_versions, num = bam_numbers, hla_gene = hla_genes)

rule hla_imputation_all:
    input:
        bamlist = expand("results/hla/imputation/bamlists/bamlist{num}.txt", num = bam_numbers),
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_v{IPD_IMGT_version}/HLA{gene}fullallelesfilledin.RData", gene = hla_genes, IPD_IMGT_version = IPD_IMGT_versions),
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
        beagle_phased_GAMCC = expand("results/phasing/HLA_GAMCC_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.GAMCC.tsv", gene = HLA_GENES, sample = samples_fv_gm)
        # vcf = expand("results/phasing/HLA_{study}_BEAGLE/unphased.{study}.chr6.vcf.gz", study = studies),
        # phased_vcf = expand("results/phasing/HLA_{study}_BEAGLE/phased.{study}.chr6.vcf.gz", study = studies),
        # oneKG_phase_df = expand("results/phasing/phased_dfs/oneKG_{vcf_version}-{filter}-{gene}.tsv", gene = HLA_GENES, filter = filters, vcf_version = vcf_versions),
        #GAMCC_phase_df = expand("results/phasing/phased_dfs/GAMCC-{gene}.tsv", gene = HLA_GENES),

        # beagle_phased_oneKG = expand("results/phasing/HLA_1KG_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.1KG.tsv", gene = HLA_GENES, sample = samples_oneKG),
        # beagle_phased_GAMCC = expand("results/phasing/HLA_GAMCC_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.GAMCC.tsv", gene = HLA_GENES, sample = samples_fv_gm),
        # concordance_df = expand("results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}.tsv", filter = filters, vcf_version = vcf_versions[0])

rule hla_ref_panel_all:
    input:
#        fv_vcf = expand("results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf.gz", chr = chromosome),
#        mGen_chunk_vcfs = [mGen_chunk_vcf_lst],
#        oneKG_gamcc = expand("results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz", chr = chromosome)
        oneKG_gamcc = expand("results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz", chr = [6]),
        oneKG_gamcc_hla = "results/hla_ref_panel/oneKG_mGenv1/merged/hla_only/oneKG_GAMCC.hla.vcf.gz"

rule hla_imputation_wip_all:
    input:
        imputed_method_v3390 = expand("results/hla/imputation/QUILT_HLA_result_method/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv)
        #ref_panel_optimal = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, id = samples_fv),
        #imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv)

        # ref_panel_optimal = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, id = samples_fv),
        # imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),

        # merged_ref_hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
        # merged_ref_legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
        # merged_ref_sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples",
        # output_db = expand('/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{hla_gene}.ssv', hla_gene = hla_genes, IPD_IMGT_version = IPD_IMGT_versions),
        # bamlist = expand("results/hla/imputation/bamlists_fv/bamlist{num}.txt", num = bam_numbers),

        # ref_panel_db = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
        # ref_panel_merged_ref_v3390 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, id = samples_fv),
        # ref_panel_merged_ref_v3570 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref_v3570/no_{id}/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, id = samples_fv),
        # ref_panel_method_v3390 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
        # ref_panel_method_v3570 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method_v3570/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, id = samples_fv),
        # ref_panel_optimal = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/hla{hla_gene}haptypes.RData", hla_gene = hla_genes, id = samples_fv),

        # imputed_db = expand("results/hla/imputation/QUILT_HLA_result_db/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers),
        # imputed_merged_ref_v3390 = expand("results/hla/imputation/QUILT_HLA_result_merged_ref/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        # imputed_merged_ref_v3570 = expand("results/hla/imputation/QUILT_HLA_result_merged_ref_v3570/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        # imputed_method_v3390 = expand("results/hla/imputation/QUILT_HLA_result_method/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        # imputed_method_v3570 = expand("results/hla/imputation/QUILT_HLA_result_method_v3570/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        # imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),

rule hla_imputation_method_all:
    input:
        db = expand('/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{gene}.ssv', gene = HLA_GENES_ALL, IPD_IMGT_version = IPD_IMGT_versions),
        db_filtered = expand('/well/band/users/rbx225/recyclable_files/hla_reference_files/v{IPD_IMGT_version}_{panel}_only/{gene}.ssv', gene = HLA_GENES_ALL, panel = panels, IPD_IMGT_version = IPD_IMGT_versions),

        # ref_panel_method_v3390 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData", gene = hla_genes),

        # imputed_method_v3390 = expand("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.combined.all.txt", gene = hla_genes, id = samples_fv)

