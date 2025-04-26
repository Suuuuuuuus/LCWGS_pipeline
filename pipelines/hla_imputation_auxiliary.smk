configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import io
import os
import re
import json
import pandas as pd
import numpy as np
import math
import pyreadr
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_test/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']

extract_dir = "results/hla/imputation/QUILT_HLA_result_method/"
rule extract_QUILT_alignments:
    input:
        RData = f"{extract_dir}{{id}}/{{gene}}/quilt.RDataput.hla{{gene}}.RData"
    output:
        extracted_RData = f"{extract_dir}{{id}}/{{gene}}/extracted.hla{{gene}}.RData"
    localrule: True
    script:
        "/well/band/users/rbx225/software/QUILT_test/QUILT/Python/extract_QUILT_alignments.R"

score_diff_in_alignment_genes_ary = [0, 8]
n_mismatches_ary = [3, 5]
weight_ary = ['T', 'F']

rule determine_imputation_hyperparameters:
    input:
        RData = f"{extract_dir}{{id}}/{{gene}}/extracted.hla{{gene}}.RData",
        bam = "data/bams/{id}.bam",
        matrix = "results/hla/imputation/WFA_alignments/v3390/{id}/{gene}/AS_matrix.ssv"
    output:
        imputed_all = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.combined.all.txt",
        imputed_top = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.combined.topresult.txt",
        quilt_all = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.onlystates.all.txt",
        quilt_top = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.onlystates.topresult.txt",
        read_all = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.onlyreads.all.txt",
        read_top = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/quilt.hla.output.onlyreads.topresult.txt"
    params:
        outputdir = "results/hla/imputation/QUILT_HLA_result_method_determine_optimal/{score}_{n_mismatches}_{weight}/{id}/{gene}/"
    run:
        os.makedirs(params.outputdir, exist_ok=True)

        d = wildcards.score
        n = wildcards.n_mismatches
        w = wildcards.weight
        if w == 'T':
            w = True 
        else:
            w = False 
        g = wildcards.gene
        s = wildcards.id

        bam = input.bam
        db_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_merged_only/'
        as_dir = 'results/hla/imputation/WFA_alignments/v3390/'
        hla_gene_information = pd.read_csv(HLA_GENE_INFORMATION_FILE, sep = ' ')

        r1, _, mate, pair = hla_aligner(g, bam, db_dir, as_dir, hla_gene_information,
                                    n_mismatches = int(n), score_diff_in_alignment_genes = int(d))

        qmat = pyreadr.read_r(input.RData)['quilt']

        a1 = set(pair.columns)
        a2 = set(qmat.columns)

        cols = sorted(list(a1.intersection(a2)))
        pair = pair.loc[cols, cols]
        qmat = qmat.loc[cols, cols]

        bestalleles = get_best_alleles(np.log(qmat))
        bestalleles.columns = ['bestallele1', 'bestallele2', 'post_prob', 'sums']
        bestalleles['sample_number'] = 1
        bestalleles['sample_name'] = s
        bestalleles = bestalleles[['sample_number', 'sample_name', 'bestallele1', 'bestallele2', 'post_prob', 'sums']]
        bestalleles.to_csv(output.quilt_all, index = False, header = True, sep = '\t')
        bestalleles.iloc[[0], :].to_csv(output.quilt_top, index = False, header = True, sep = '\t')

        pair = pair - pair.max().max()
        mask = np.tril(np.ones(pair.shape), k=-1).astype(bool)
        pair = pair.where(~mask, other=0)
        pair = np.exp(pair)
        pair = pair/pair.sum().sum()

        bestalleles = get_best_alleles(np.log(pair))
        bestalleles.columns = ['bestallele1', 'bestallele2', 'post_prob', 'sums']
        bestalleles['sample_number'] = 1
        bestalleles['sample_name'] = s
        bestalleles = bestalleles[['sample_number', 'sample_name', 'bestallele1', 'bestallele2', 'post_prob', 'sums']]
        bestalleles.to_csv(output.read_all, index = False, header = True, sep = '\t')
        bestalleles.iloc[[0], :].to_csv(output.read_top, index = False, header = True, sep = '\t')

        if w:
            w1 = pair.max().max()
            w2 = qmat.max().max()
            w1 = w1/(w1+w2)
            w2 = 1-w1
            
            combined = (pair**w1)*(qmat**w2)
        else:
            combined = pair*qmat
            
        combined = combined/combined.sum().sum()

        bestalleles = get_best_alleles(np.log(combined))
        bestalleles.columns = ['bestallele1', 'bestallele2', 'post_prob', 'sums']
        bestalleles['sample_number'] = 1
        bestalleles['sample_name'] = s
        bestalleles = bestalleles[['sample_number', 'sample_name', 'bestallele1', 'bestallele2', 'post_prob', 'sums']]

        bestalleles.to_csv(output.imputed_all, index = False, header = True, sep = '\t')
        bestalleles.iloc[[0], :].to_csv(output.imputed_top, index = False, header = True, sep = '\t')