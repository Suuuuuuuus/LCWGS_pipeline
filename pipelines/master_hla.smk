include: "hla.smk"
#include: "alignment.smk"
#include: "post_hla.smk"
#include: "hla_ref_panel.smk"
include: "hla_imputation_wip.smk"
include: "hla_imputation_prep.smk"
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
samples_oneKG = read_tsv_as_lst("/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/samples_to_phase.tsv")
chromosome = [i for i in range(1,23)]
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
studies = ['1KG', 'GAMCC']
filters = ['strict', 'loose']
vcf_versions = ['30x', 'phase3_b38']

rule hla_all:
    input:
        called = expand("results/hla/call/{id}/hla/R1_bestguess_G.txt", id = samples_fv)

rule hla_imputation_prep_all:
    input:
        reads1 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv", hla_gene = HLA_GENES, id = samples_fv),
        reads2 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv", hla_gene = HLA_GENES, id = samples_fv),
        mate_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv", hla_gene = HLA_GENES, id = samples_fv),
        pair_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv", hla_gene = HLA_GENES, id = samples_fv)
        

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
        #vcf = "results/phasing/HLA_1KG_BEAGLE/unphased.1KG.chr6.vcf.gz",
        #phased_vcf = "results/phasing/HLA_1KG_BEAGLE/phased.1KG.chr6.vcf.gz",
        # oneKG_html = expand("results/phasing/html/oneKG-{filter}-{gene}.html", gene = HLA_GENES, filter = filters), 
        oneKG_phase_df = expand("results/phasing/phased_dfs/oneKG_{vcf_version}-{filter}-{gene}.tsv", gene = HLA_GENES, filter = filters, vcf_version = vcf_versions),
        # GAMCC_html = expand("results/phasing/html/GAMCC-{gene}.html", gene = HLA_GENES), 
        #GAMCC_phase_df = expand("results/phasing/phased_dfs/GAMCC-{gene}.tsv", gene = HLA_GENES),

        beagle_phased = expand("results/phasing/HLA_1KG_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.1KG.tsv", gene = HLA_GENES, sample = samples_oneKG),
        concordance_df = expand("results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}.tsv", filter = filters, vcf_version = vcf_versions[0])

rule hla_ref_panel_all:
    input:
        # hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        # legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        # samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples",
        # RData = "results/hla_tests/quilt.hrc.hla.all.haplotypes.RData",
        # bam_all = "results/hla_tests/bamlist.txt",
        fv_vcf = expand("results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf.gz", chr = chromosome),
        mGen_chunk_vcfs = [mGen_chunk_vcf_lst],
        oneKG_gamcc = expand("results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz", chr = chromosome)

rule hla_imputation_wip_all:
    input:
        output_db = expand('/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv', hla_gene = hla_genes),
        bamlist = expand("results/hla/imputation/bamlists_fv/bamlist{num}.txt", num = bam_numbers),
#        merged_ref_hap = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.hap.gz",
#        merged_ref_legend = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.legend.gz",
#        merged_ref_sample = "results/hla_ref_panel/oneKG_mGenv1/hla_prepare_ref/oneKG_GAMCC.chr6.samples",

        ref_panel_db = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
#        ref_panel_merged_ref = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_merged_ref/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, id = samples_fv),
        ref_panel_method = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes),
#        ref_panel_optimal = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_optimal/no_{id}/HLA{hla_gene}fullallelesfilledin.RData", hla_gene = hla_genes, id = samples_fv),

        imputed_db = expand("results/hla/imputation/QUILT_HLA_result_db/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers),
#        imputed_merged_ref = expand("results/hla/imputation/QUILT_HLA_result_merged_ref/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),
        imputed_method = expand("results/hla/imputation/QUILT_HLA_result_method/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, num = bam_numbers),
#        imputed_optimal = expand("results/hla/imputation/QUILT_HLA_result_optimal/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = hla_genes, id = samples_fv),

