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

samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
samples_fv_gam = read_tsv_as_lst('data/sample_tsvs/fv_gam_names.tsv')
chromosome = [i for i in range(1,23)]
bam_batches = config['bam_batch']
bam_numbers = [str(i) for i in range(1, int(bam_batches) + 1)]
hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
v3570_db_alleles = [5093, 6106, 5657, 869, 624]

rule hla_imputation_prepare_align_per_allele:
    input:
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        read1_ary = temp("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary1.npy"),
        read2_ary = temp("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary2.npy")
    resources:
        mem = '40G'
    threads: 4
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv"
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        db = pd.read_csv(input.db_file, sep = ' ')
        gene = wildcards.hla_gene
        reads_apart_max = 10000
        bam = input.bam

        reads1 = get_chr6_reads(gene, bam, hla_gene_information, reads_apart_max)
        reads2 = get_hla_reads(gene, bam, reads_apart_max)

        if reads1.empty:
            reads1 = reads2.iloc[:2, :] if not reads2.empty else pd.DataFrame()
        elif reads2.empty:
            reads2 = reads1.iloc[:2, :]
        else:
            pass

        reads1['rev_seq'] = reads1['sequence'].apply(reverse_complement)
        reads1['rev_bq'] = reads1['base_quality'].apply(lambda bq: bq[::-1])

        j = wildcards.idx
        a = db.columns[j]
        res = per_allele(j, a, db, reads1)
        scores_ary = np.maximum(res[1], res[2])

        np.set_printoptions(precision=6, suppress=True)
        np.save(output.read1_ary, scores_ary)

        reads2['rev_seq'] = reads2['sequence'].apply(reverse_complement)
        reads2['rev_bq'] = reads2['base_quality'].apply(lambda bq: bq[::-1])

        res = per_allele(j, a, db, reads2)
        scores_ary = np.maximum(res[1], res[2])
        
        np.set_printoptions(precision=6, suppress=True)
        np.save(output.read2_ary, scores_ary)

rule hla_imputation_prepare_per_sample:
    input:
        read1_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary1.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        read2_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary2.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        reads1 = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv",
        reads2 = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv",
        mate_matrix = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv",
        pair_matrix = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv"
    resources:
        mem = '60G'
    threads: 4
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv"
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        db = pd.read_csv(input.db_file, sep = ' ')
        gene = wildcards.hla_gene
        reads_apart_max = 10000
        bam = input.bam

        reads1 = get_chr6_reads(gene, bam, hla_gene_information, reads_apart_max)
        reads2 = get_hla_reads(gene, bam, reads_apart_max)

        scores_mat1 = [np.hstack([np.load(f'results/hla/imputation/ref_panel/QUILT_bam_aligner/{wildcards.hla_gene}/{wildcards.id}/{wildcards.id}.{idx}.ary1.npy')[:,np.newaxis]]) for idx in range(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)])]
        likelihood_mat1 = np.exp(scores_mat1)/np.sum(np.exp(scores_mat1), axis = 1, keepdims = True)
        likemat1 = np.log(likelihood_mat1)

        scores_mat2 = [np.hstack([np.load(f'results/hla/imputation/ref_panel/QUILT_bam_aligner/{wildcards.hla_gene}/{wildcards.id}/{wildcards.id}.{idx}.ary2.npy')[:,np.newaxis]]) for idx in range(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)])]
        likelihood_mat2 = np.exp(scores_mat2)/np.sum(np.exp(scores_mat2), axis = 1, keepdims = True)
        likemat2 = np.log(likelihood_mat2)

        rl = reads1['sequence'].str.len().mode().values[0]
        n_mismatches = 5
        assumed_bq = 0.001
        min_valid_prob = np.log(math.comb(rl, n_mismatches)) + n_mismatches*np.log(assumed_bq) + (rl - n_mismatches)*np.log(1 - assumed_bq)

        valid_indices1 = np.any(likemat1 >= min_valid_prob, axis=1)
        valid_indices2 = np.any(likemat2 >= min_valid_prob, axis=1)
        likemat1, reads1 = likemat1[valid_indices1], reads1[valid_indices1]
        likemat2, reads2 = likemat2[valid_indices2], reads2[valid_indices2]
        
        likemat_all = np.vstack((likemat1, likemat2))

        id1, id2 = reads1.iloc[:, 0].to_numpy(), reads2.iloc[:, 0].to_numpy()
        
        readind = (reads1.iloc[:, 1].astype(int) // 64) % 4
        readind2 = (reads2.iloc[:, 1].astype(int) // 64) % 4
        mate_indicator = np.concatenate((readind, readind2))

        ids_all = np.concatenate((id1, id2))
        unique_ids = np.unique(ids_all)
        likemat_mate = np.zeros((len(unique_ids), likemat_all.shape[1]))

        for i, uid in enumerate(unique_ids):
            t1 = likemat_all[ids_all == uid, :]
            t2 = mate_indicator[ids_all == uid]
            if len(t2) > 0:
                likemat_mate[i, :] = np.sum(t1[t2 > 0], axis=0)

        valid_mask = likemat_mate.max(axis=1) >= min_valid_prob
        likemat_mate = likemat_mate[valid_mask]
        likemat_norm = 0.5 * np.exp(likemat_mate - likemat_mate.max(axis=1, keepdims=True)) + 1e-100

        likemat_paired = likemat_norm.T @ likemat_norm
        likemat_paired = pd.DataFrame(likemat_paired, index=db.columns, columns=db.columns)
        likemat_mate = pd.DataFrame(likemat_mate, index = unique_ids[valid_mask], columns=db.columns)

        reads1.to_csv(output.reads1, header = False, index = False)
        reads2.to_csv(output.reads2, header = False, index = False)
        pd.set_option('display.float_format', '{:.6e}'.format)
        likemat_mate.to_csv(output.mate_matrix, index=True, header=True, sep = ' ')
        likemat_paired.to_csv(output.pair_matrix, index=True, header=True, sep = ' ')

            