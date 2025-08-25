import io
import os
import sys
import csv
import gzip
import time
import json
import secrets
import copy
import pickle
import multiprocessing
import subprocess
import resource
import itertools
from itertools import combinations_with_replacement
import collections
import sqlite3
import random

import pyranges as pr
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import scipy as sp
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

from scipy.stats import nbinom
from scipy.stats import geom, beta
from scipy.special import logsumexp

sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
import lcwgsus
from lcwgsus.variables import *

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

def simulate_coverage(L, Ecov, Dcov, binsize=1000, method="negative_binomial"):
    Ecov = np.where(Ecov * binsize == 0, 1e-3, Ecov * binsize)
    Dcov = np.where(Dcov * binsize == 0, 1e-2, Dcov * binsize)

    if method == "negative_binomial":
        size = (Ecov**2) / (Dcov - Ecov)
        prob = Ecov / Dcov
        return nbinom.rvs(size, prob, size=L)
    elif method == "normal":
        return np.maximum(norm.rvs(loc=Ecov, scale=np.sqrt(Dcov) + 0.05, size=S), 0)
    elif method == "poisson":
        return poisson.rvs(mu=Ecov, size=L)
    else:
        raise ValueError("Select from: 1.normal; 2.negative_binomial; 3.poisson.")

def simulate_coverotron(model, N, L, mean_coverage, sd_coverage, nb_var, binsize = 1000, model_sim = 'negative_binomial', ploidy = 2):
    cumulative = np.cumsum(np.concatenate([[0], model.freqs]))

    means = np.random.normal(loc=mean_coverage, scale=sd_coverage, size=N)
    variances = nb_var * means

    random_vals = np.random.rand(ploidy * N)
    hap_indices = [np.where(cumulative < val)[0][-1] for val in random_vals]
    haplotypes = np.array(hap_indices).reshape(N, ploidy)

    profiles = np.full((L, N), np.nan)
    coverage = np.full((L, N), np.nan)
    training = np.full((L, N), np.nan)

    for i in range(N):
        training[:, i] = simulate_coverage(L, means[i] * ploidy, variances[i] * ploidy, binsize, method=model_sim)

        profile_1 = model.haps[haplotypes[i, 0]]
        profile_2 = model.haps[haplotypes[i, 1]]
        profiles[:, i] = profile_1 + profile_2

        Ecov = means[i] * profiles[:, i]
        Dcov = variances[i] * profiles[:, i]
        coverage[:, i] = simulate_coverage(L, Ecov, Dcov, binsize, method=model_sim)

    haplotypes = haplotypes.astype(int)
    haplotypes = [tuple(sorted(i)) for i in haplotypes]
    return training, coverage, haplotypes

def deresolute_windows(df, window_size, normalise = False):
    len_window = int(df.iloc[1,0] - df.iloc[0,0])
    len_region = int(len(df)*len_window)
    if window_size % len_window != 0 or (len_region % window_size) != 0:
        raise ValueError(f"Window size must be a multiple of {len_window} and divide {len_region}.")
    multiple = window_size/len_window
        
    df['window_start'] = (df.iloc[:, 0] // window_size) * window_size
    new_pos = df['window_start'].unique()
    agg_df = df.groupby('window_start').sum().reset_index()
    agg_df['position'] = new_pos
    agg_df = agg_df.drop(columns=['window_start'])
    
    if normalise:
        for s in agg_df.columns[1:]:
            agg_df[s] = agg_df[s]/window_size
    return agg_df

def delineate_region(start, length, binsize = 1000):
    s = start - length
    e = start + 2*length
    
    s = int(s/binsize)*binsize
    e = (int(e/binsize) + 1)*binsize
    return s,e

def delineate_region2(eichler_full, c, start, length, svtype, binsize = 1000, extension = 3):
    if svtype == 'INS':
        e = start + length
        s = start - length
    else:
        e = start + length
        s = start
        
    s,e = find_overlapping_svs(eichler_full, c, s, e)
    
    s = int(s/binsize)*binsize - extension*binsize
    e = (int(e/binsize) + 1)*binsize + extension*binsize
    return s,e

def delineate_region3(start, length, svtype, binsize = 1000, extension = 3):
    if svtype == 'INS':
        e = start + length
        s = start - length
    else:
        e = start + length
        s = start
    
    s = int(s/binsize)*binsize - extension*binsize
    e = (int(e/binsize) + 1)*binsize + extension*binsize
    return s,e

def read_eichler(f, length_filter = 0, af_filter = 0):
    df = pd.read_csv(f, sep = '\t', compression = 'gzip')
    main_chrs = [f'chr{i}' for i in range(1,23)] + ['X', 'Y']
    df = df[df['#CHROM'].isin(main_chrs)]
    df['PG_AFR_AF'] =  df['PG_INFO_AFR'].str.split(';').str.get(0).str.split('=').str.get(1).astype(float)
    df = df.sort_values(by = 'POP_AFR_AF', ascending = False)
    df = df[(df['SVLEN'] >= length_filter)]
#     df = df[(df['SVLEN'] >= length_filter) & ((df['POP_AFR_AF'] >= af_filter) | (df['PG_AFR_AF'] >= af_filter))]
    df = df[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN','PG_AFR_AF', 'POP_AFR_AF']]
    df = df.sort_values(by = df.columns[:3].tolist()).reset_index(drop = True)
    return df

def parse_eichler_full(f = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz'):
    eichler_full = read_eichler(f)
    indices = np.where(eichler_full['SVTYPE'] == 'INS')[0]
    eichler_full.loc[indices, 'END'] = eichler_full.loc[indices, 'POS'] + eichler_full.loc[indices, 'SVLEN']
    eichler_full.loc[indices, 'POS'] = eichler_full.loc[indices, 'POS'] - eichler_full.loc[indices, 'SVLEN']
    cols = ['Chromosome', 'Start', 'End']
    eichler_full.columns = cols + eichler_full.columns[3:].tolist()
    eichler_full = eichler_full.sort_values(by = cols).reset_index(drop = True)
    return eichler_full

def find_eichler_region(ix, eichler_full, eichler_manifest = 'results/nonahore/eichler/manifest.5K.5percent.tsv'):
    manifest = pd.read_csv(eichler_manifest, sep = '\t')
    chrom, sv_start, L, svtype = manifest.loc[ix, ['#CHROM', 'POS', 'SVLEN', 'SVTYPE']]
    chromosome = int(chrom.replace('chr', ''))
    start, end = delineate_region2(eichler_full, chrom, sv_start, L, svtype)
    return start, end

def parse_nonahore_result1(bed = '/well/band/users/rbx225/GAMCC/data/bedgraph/GRCh38.autosomes.bed', 
                          eichler_manifest = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
                         eichler_vcf = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel_alt.vcf.gz'):
    df = pd.read_csv(eichler_manifest, sep = '\t', compression = 'gzip')
    main_chrs = [f'chr{i}' for i in range(1,23)] + ['X', 'Y']
    df = df[df['#CHROM'].isin(main_chrs)]
    df['PG_AFR_AF'] =  df['PG_INFO_AFR'].str.split(';').str.get(0).str.split('=').str.get(1).astype(float)
    df = df.sort_values(by = 'POP_AFR_AF', ascending = False)
    df = df[(df['SVLEN'] > 5000) & ((df['POP_AFR_AF'] >= 0.05) | (df['PG_AFR_AF'] >= 0.05))]
    df = df[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN','ID']].reset_index(drop = True)

    x = lcwgsus.read_vcf(eichler_vcf)
    x = x[x['ID'].isin(df['ID'])].reset_index(drop = True)
    df = pd.merge(df, x[['ID', 'ref', 'alt', 'INFO']], on = ['ID'])
    df['TRF'] = df['INFO'].apply(lambda x: 1 if 'REF_TRF' in x else 0)

    n_hap_ary = []
    true_info_ary = []
    concordance_ary = []
    est_freq_ary = []
    mq_ary = []
    bad_alignment_region_ary = []
    true_hap_idx_ary = []
    second_info_ary = []
    max_info_ary = []
    max_maf_ary = []

    for ix in range(len(df)):
        infile = f'results/nonahore/eichler/nonahore1/region{ix}/results.pickle'
        if os.path.exists(infile):
            data = read_pickle(infile)
            n_hap = len(data['haps'])
            n_hap_ary.append(n_hap)
            est_freq_ary.append(data['freq'])
            concordance_ary.append(data['concordance'])
            mq_ary.append(data['mean_mq'])
            bad_alignment_region_ary.append(check_low_mq(data['mq']))

            true_hap_idx = data['true_hap_idx']
            if true_hap_idx == -9:
                true_hap_idx_ary.append(-9)
                true_info_ary.append(0)
            else:
                true_hap_idx_ary.append(true_hap_idx)
                true_info_ary.append(data['info'][true_hap_idx])

            if n_hap == 1:
                max_info_ary.append(0)
                max_maf_ary.append(0)
                second_info_ary.append(0)
            else:
                max_info_idx = np.array(data['info'][1:]).argmax() + 1
                max_info_ary.append(data['info'][max_info_idx])
                max_maf_ary.append(data['freq'][max_info_idx])
                second_info_ary.append(data['info'][1])
        else:
            n_hap_ary.append(0)
            true_info_ary.append(0)
            est_freq_ary.append(0)     
            concordance_ary.append(0)
            mq_ary.append(0)
            bad_alignment_region_ary.append(0)
            true_hap_idx_ary.append(-9)
            max_info_ary.append(0)
            max_maf_ary.append(0)
            second_info_ary.append(0)

    df['n_hap'] = n_hap_ary
    df['freq'] = est_freq_ary
    df['concordance'] = concordance_ary
    df['mean_mq'] = mq_ary
    df['contain_bad_alignment'] = bad_alignment_region_ary

    df['true_hap_idx'] = true_hap_idx_ary
    df['true_info'] = true_info_ary
    df['max_info'] = max_info_ary
    df['max_maf'] = max_maf_ary
    df['second_info'] = second_info_ary

#     bed = pd.read_csv(bed, sep = '\t', header = None)
#     bed.columns = ['chr', 'start', 'end']

    df1 = df.copy()
#     df1['is_endchr'] = False
#     df1 = df1.apply(check_endchr, axis = 1)
    return df1

def parse_nonahore_result2(bed = '/well/band/users/rbx225/GAMCC/data/bedgraph/GRCh38.autosomes.bed', 
                          eichler_manifest = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel.tsv.gz',
                         eichler_vcf = '/well/band/users/rbx225/recyclable_files/eichler_sv/variants_freeze4_sv_insdel_alt.vcf.gz'):
    df = pd.read_csv(eichler_manifest, sep = '\t', compression = 'gzip')
    main_chrs = [f'chr{i}' for i in range(1,23)] + ['X', 'Y']
    df = df[df['#CHROM'].isin(main_chrs)]
    df['PG_AFR_AF'] =  df['PG_INFO_AFR'].str.split(';').str.get(0).str.split('=').str.get(1).astype(float)
    df = df.sort_values(by = 'POP_AFR_AF', ascending = False)
    df = df[(df['SVLEN'] > 5000) & ((df['POP_AFR_AF'] >= 0.05) | (df['PG_AFR_AF'] >= 0.05))]
    df = df[['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN','ID']].reset_index(drop = True)

    x = lcwgsus.read_vcf(eichler_vcf)
    x = x[x['ID'].isin(df['ID'])].reset_index(drop = True)
    df = pd.merge(df, x[['ID', 'ref', 'alt', 'INFO']], on = ['ID'])
    df['TRF'] = df['INFO'].apply(lambda x: 1 if 'REF_TRF' in x else 0)

    n_hap_ary = []
    true_info_ary = []
    concordance_ary = []
    est_freq_ary = []
    mq_ary = []
    bad_alignment_region_ary = []
    true_hap_idx_ary = []
    second_info_ary = []
    max_info_ary = []
    max_maf_ary = []
    to_call_ary = []
    callable_windows_ary = []

    for ix in range(len(df)):
        infile = f'results/nonahore/eichler/nonahore2/region{ix}/results.pickle'
        if os.path.exists(infile):
            data = read_pickle(infile)
            to_call = data['to_call']
            callable_windows_ary.append(data['callable_windows'])
        else:
            data = {}
            to_call = False
            callable_windows_ary.append(0)
        
        to_call_ary.append(to_call)
            
        if to_call:
            n_hap = len(data['haps'])
            n_hap_ary.append(n_hap)
            est_freq_ary.append(data['freq'])
            concordance_ary.append(data['concordance'])
            mq_ary.append(data['mean_mq'])
            bad_alignment_region_ary.append(check_low_mq(data['mq']))

            true_hap_idx = data['true_hap_idx']
            if true_hap_idx == -9:
                true_hap_idx_ary.append(-9)
                true_info_ary.append(0)
            else:
                true_hap_idx_ary.append(true_hap_idx)
                true_info_ary.append(data['info'][true_hap_idx])

            if n_hap == 1:
                max_info_ary.append(0)
                max_maf_ary.append(0)
                second_info_ary.append(0)
            else:
                max_info_idx = np.array(data['info'][1:]).argmax() + 1
                max_info_ary.append(data['info'][max_info_idx])
                max_maf_ary.append(data['freq'][max_info_idx])
                second_info_ary.append(data['info'][1])
        else:
            n_hap_ary.append(0)
            true_info_ary.append(0)
            est_freq_ary.append(0)     
            concordance_ary.append(0)
            mq_ary.append(0)
            bad_alignment_region_ary.append(0)
            true_hap_idx_ary.append(-9)
            max_info_ary.append(0)
            max_maf_ary.append(0)
            second_info_ary.append(0)

    df['n_hap'] = n_hap_ary
    df['freq'] = est_freq_ary
    df['concordance'] = concordance_ary
    df['mean_mq'] = mq_ary
    df['contain_bad_alignment'] = bad_alignment_region_ary

    df['true_hap_idx'] = true_hap_idx_ary
    df['true_info'] = true_info_ary
    df['max_info'] = max_info_ary
    df['max_maf'] = max_maf_ary
    df['second_info'] = second_info_ary
    
    df['to_call'] = to_call_ary
    df['callable_windows'] = callable_windows_ary

#     bed = pd.read_csv(bed, sep = '\t', header = None)
#     bed.columns = ['chr', 'start', 'end']

    df1 = df.copy()
#     df1['is_endchr'] = False
#     df1 = df1.apply(check_endchr, axis = 1)
    return df1

def read_tabix(vcf, c, regstart, regend):
    command = f"tabix {vcf} chr{c}:{regstart}-{regend}"
    seqs = subprocess.run(command, shell = True, capture_output = True, text = True).stdout[:-1].split('\n')
    if seqs == ['']:
        return pd.DataFrame()

    seqs = [i.split('\t') for i in seqs if '##' not in i]
    seqs = pd.DataFrame(seqs).iloc[:, :5]
    seqs.columns = ['chr', 'pos', 'id', 'ref', 'alt']
    return seqs

def read_bed_tabix(vcf, c, regstart, regend):
    command = f"tabix {vcf} chr{c}:{regstart}-{regend}"
    seqs = subprocess.run(command, shell = True, capture_output = True, text = True).stdout[:-1].split('\n')
    if seqs == ['']:
        return pd.DataFrame()

    seqs = pd.DataFrame([i.split('\t') for i in seqs if '##' not in i], 
                        columns=["chr", "start", "end", "value"]).astype({"start": int, "end": int, "value": float})
    
    if seqs.loc[0, 'start'] < regstart:
        seqs.loc[0, 'start'] = regstart
    elif seqs.loc[0, 'start'] > regstart:
        new_row = {'chr': c, 'start': regstart, 'end': seqs.loc[0, 'start'], 'value': 0}
        seqs = pd.concat([pd.DataFrame([new_row]), seqs], ignore_index=True).reset_index(drop = True)
    else:
        pass

    if seqs.loc[len(seqs)-1, 'end'] > regend:
        seqs.loc[len(seqs)-1, 'end'] = regend
    elif seqs.loc[len(seqs)-1, 'end'] < regend:
        new_row = {'chr': c, 'start': seqs.loc[len(seqs)-1, 'end'], 'end': regend, 'value': 0}
        seqs = pd.concat([seqs, pd.DataFrame([new_row])], ignore_index=True).reset_index(drop = True)
    else:
        pass
    
    filled = []
    prev_end = regstart
    for _, row in seqs.iterrows():
        if row["start"] > prev_end:
            filled.append([prev_end, row["start"], 0.0])
        filled.append([row["start"], row["end"], row["value"]])
        prev_end = row["end"]

    if prev_end < regend:
        filled.append([prev_end, regend, 0.0])
        
    filled = pd.DataFrame(filled, columns = seqs.columns[1:])
    filled['chr'] = c
    filled = filled[['chr'] + seqs.columns[1:].tolist()]
    return filled

def read_bam(bam, c, regstart, regend):
    command = f"samtools view {bam} chr{c}:{regstart}-{regend}"
    reads = subprocess.run(command, shell = True, capture_output = True, text = True).stdout[:-1].split('\n')
    if reads == ['']:
        return pd.DataFrame()

    reads = [i.split('\t') for i in reads if '##' not in i]
    reads = pd.DataFrame(reads)
    reads.columns = [
        "ID", "flag", "chr", "pos", "map_quality", "CIGAR", "chr_alt", "pos_alt", "insert_size", "sequence", "base_quality"
    ] + [f'col{i}' for i in range(11, reads.shape[1])]

    reads['pos'] = reads['pos'].astype(int)
    reads['pos_alt'] = reads['pos_alt'].astype(int)
    return reads

def read_pickle(infile):
    with open(infile, 'rb') as f:
        data = pickle.load(f)
    return data

def read_coverage_data(file_path, sep = ','):    
    df = pd.read_csv(file_path, sep = sep)
    bins = list(zip(df['position'], df['position'] + df['size'], df['N']))
    samples = [col for col in df.columns if col.endswith(":coverage") and col != "total:coverage"]
    coverage = df.loc[:,samples].replace("NA", 0).astype(float).to_numpy()
    samples = [s.split(':')[0] for s in samples]
    return samples, bins, coverage

def check_endchr(r, cutoff = 1e6):
    c = r['#CHROM']
    s = r['POS']
    endchr = bed[bed['chr'] == c]['end'].values[0]

    if s > cutoff and s < endchr - cutoff:
        r['is_endchr'] = False
    else:
        r['is_endchr'] = True
    return r

def calculate_mappability_per_bin(df, start, end, binsize=1000):
    bins = np.arange(start, end, binsize)
    n_bins = len(bins)
    out = pd.DataFrame({
        "start": bins,
        "end": bins + binsize,
        "total_score": 0.0,
        "total_len": 0
    })

    for _, row in df.iterrows():
        istart, iend, val = row["start"], row["end"], row["value"]
        b0 = (istart - start) // binsize
        b1 = min((iend - start) // binsize, n_bins - 1)
        for b in range(int(b0), int(b1) + 1):
            bstart, bend = out.loc[b, ["start", "end"]]
            overlap = max(0, min(iend, bend) - max(istart, bstart))
            if overlap > 0:
                out.at[b, "total_score"] += val * overlap
                out.at[b, "total_len"]   += overlap

    out["mappability"] = out["total_score"] / out["total_len"]
    out = out[["start","mappability"]]
    out.columns = ['position', 'mappability']
    return out

def assess_bins(df, plausible_boundaries, map_t = 0.9, mq_t = 20, percent_windows_t = 0.5):
    start_idx, end_idx = plausible_boundaries
    metrics = df.copy()
    indices = np.where((metrics['mappability'] >= map_t) & (metrics['total:mean_mq'] >= mq_t))[0]
    metrics['valid'] = 0
    metrics.loc[indices, 'valid'] = 1
    
    tmp = metrics.loc[start_idx:end_idx,:]
    percent_windows = tmp['valid'].sum()/len(tmp)
    if percent_windows <= percent_windows_t:
        to_call = False
    else:
        to_call = True
    include_bins = np.where(metrics['valid'] == 1)[0]
    return to_call, include_bins

def extract_target_cov(df, start, end):
    df = df[(df['position'] >= start) & (df['position'] < end)].reset_index(drop = True)
    samples = [col for col in df.columns if col.endswith(":coverage") and col != "total:coverage"]
    coverage = df.loc[:,samples].replace("NA", 0).astype(float).to_numpy()
    samples = [s.split(':')[0] for s in samples]
    return samples, coverage

def normalise_by_flank(df, start, end, flank, side = 'both'):
    fstart = max(start-flank, df.iloc[0,0])
    fend = min(end+flank, df.iloc[-1,0])

    if side == 'both':
        criteria = (((df['position'] >= fstart) & (df['position'] < start)) | 
             ((df['position'] >= end) & (df['position'] < fend)))
    elif side == 'left':
        criteria = ((df['position'] >= fstart) & (df['position'] < start))
    elif side == 'right':
        criteria = ((df['position'] >= end) & (df['position'] < fend))
    else:
        raise ValueError("Unsupported side.")            
            
    cov = df[criteria].iloc[:,1:-1].to_numpy()
    means = np.mean(cov, axis = 0)
    variances = np.var(cov, axis = 0, ddof = 1)
    return means, variances

def normalise_by_flank2(df, start, end, flank, side = 'both'):
    fstart = max(start-flank, df.iloc[0,0])
    fend = min(end+flank, df.iloc[-1,0])
    
    left_flank = df[(df['position'] >= fstart) & (df['position'] < start)]
    right_flank = df[(df['position'] >= end) & (df['position'] < fend)]

    left_flank = left_flank[left_flank['total:coverage'] != 0] # Remove unaligned sites
    right_flank = right_flank[right_flank['total:coverage'] != 0]
    
    if side == 'both':
        cov = pd.concat([left_flank, right_flank], axis = 0)
        cov = cov.iloc[:,1:-1].to_numpy()
    elif side == 'left':
        cov = left_flank.iloc[:,1:-1].to_numpy()
    elif side == 'right':
        cov = right_flank.iloc[:,1:-1].to_numpy()
    else:
        raise ValueError("Unsupported side.")    
        
    means = np.mean(cov, axis = 0)
    variances = np.var(cov, axis = 0, ddof = 1)
    return means, variances

def load_region_files(chunks, chromosome, start, end, indir = 'results/coverage/coverotron/', flank = None):
    chromosome = str(chromosome)
    if flank is not None:
        start = start - flank
        end = end + flank
        
    starts = chunks[chromosome]['start']
    ends = chunks[chromosome]['end']
    
    dfs = []
    for s, e in zip(starts, ends):
        if not (e < start or s > end):  
            file_path = os.path.join(indir, f'chr{chromosome}.{s}.{e}.tsv.gz')
            if os.path.exists(file_path):
                dfs.append(pd.read_csv(file_path, sep='\t', compression = None))
            else:
                print(f"Warning: File {file_path} not found.")

    if dfs:
        return pd.concat(dfs, ignore_index=True)
    else:
        raise ValueError("No overlapping files found for the specified region.")

def precompute_site_lls(means, variances, coverage, max_cnv=20, mismap_proportion = 0.01, model='negative_binomial'):
    L, N = coverage.shape
    site_factor = np.ones(L)
    mismap_propn = np.full(L, mismap_proportion)

    LOG_TWO_PI = np.log(2 * np.pi)
    
    lls = []

    for cn in range(max_cnv):
        mean_term = (
            (means * cn)[np.newaxis, :] * site_factor[:, np.newaxis] +
            means[np.newaxis, :] * mismap_propn[:, np.newaxis]
        )

        if model == 'normal':
            var_term = (
                (variances * cn)[np.newaxis, :] * site_factor[:, np.newaxis] +
                variances[np.newaxis, :] * mismap_propn[:, np.newaxis]
            )
            var_term = np.maximum(var_term, 1e-10)
            residuals = coverage - mean_term
            log_likelihood = (
                -0.5 * LOG_TWO_PI
                - 0.5 * np.log(var_term)
                - 0.5 * (residuals ** 2) / var_term
            )

        elif model == 'negative_binomial':
            mu = mean_term
            var = (
                (variances * cn)[np.newaxis, :] * site_factor[:, np.newaxis] +
                variances[np.newaxis, :] * mismap_propn[:, np.newaxis]
            )
            var = np.maximum(var, mu + 1e-6)
            r = (mu ** 2) / (var - mu)
            p = r / (r + mu)
            x = np.clip(np.round(coverage), 0, None).astype(int)
            log_likelihood = nbinom.logpmf(x, r, p)

        else:
            raise ValueError("Unsupported model. Choose from: 1.normal; 2.negative_binomial.")
            
        lls.append(log_likelihood)
    return lls

class SVModel:
    def __init__(self, haps=None, freqs=None):
        self.haps = haps
        self.freqs = freqs
    
    def __repr__(self):
        return str(len(self.haps))
    
    def normalise(self):
        tmp = np.array(self.freqs)
        self.freqs = list(tmp/tmp.sum())
    
    def add(self, hap, freq):
        self.haps.append(hap)
        self.freqs = list((1-freq)*np.array(self.freqs))
        self.freqs.append(freq)
        self.normalise()
        
    def replace(self, hap, ix):
        self.haps[ix] = hap
        
def dirichlet_sampling(freqs, concentration = 10, limit = 1e-4):
    alpha = np.array(freqs) * concentration
    alpha[alpha <= limit] = limit
    rng = default_rng()
    res = rng.dirichlet(alpha)
    res[res <= limit] = limit
    return list(res)

def sample_from_freqs(freqs):
    freqs = np.array(freqs)
    cum_probs = np.cumsum(freqs)
    u = np.random.uniform(0, 1)
    return np.searchsorted(cum_probs, u)

def sample_recombinants(model, L, max_cnv = 10, check_noisy = True):
    hap1 = sample_from_freqs(model.freqs)
    hap2 = sample_from_freqs(model.freqs)
    breakpoints = np.random.choice(np.arange(1, L), size=2, replace=True)
    
    left = model.haps[hap1].copy()
    left[breakpoints[0]:] = 0
    right = model.haps[hap2].copy()
    right[:breakpoints[1]] = 0

    new_hap = left + right
    
    if np.any(new_hap >= max_cnv):
        new_hap = np.ones(L)
        
    if check_noisy:
        if not check_noisy_hap(new_hap):
            new_hap = np.ones(L)
    return new_hap

def check_noisy_hap(h, t=2):
    h = np.asarray(h)
    n = len(h)
    i = 0

    while i < n:
        if h[i] == 1:
            i += 1
            continue

        # Start of a non-1 segment
        this_sv = h[i]
        start = i
        while i < n and h[i] == this_sv:
            i += 1
        length = i - start

        if length <= t:
            return False  # Too short even though same value
    return True

def run_inverse_length_penalty(hap):
    penalty = 0
    run_length = 1
    for i in range(1, len(hap)):
        if hap[i] == hap[i - 1]:
            run_length += 1
        else:
            penalty += 1 / run_length
            run_length = 1
    penalty += 1 / run_length
    return penalty

def multi_evaluate_model(model, N, pre_computed_lls, geom_penalty, L = None, bin_size = 1000):
    n_haps = len(model.haps)
    if L is None:
        L = len(model.haps[0])
    n_diploid = int(n_haps*(n_haps + 1)/2)
    
    dip_model = generate_diploid_profiles(model)
    final_lls = np.zeros((N, n_diploid))
    
    for i in range(n_diploid):
        ll_ary = evaluate_per_hap(dip_model.haps[i], pre_computed_lls)
        final_lls[:, i] = ll_ary
        
    trans_penalty = 0
    for h in model.haps:
        trans_penalty += run_inverse_length_penalty(h)

    final_lls = final_lls + np.log(np.array(dip_model.freqs))[np.newaxis, :] 
    model_ll = logsumexp(final_lls, axis=1).sum()
    model_ll = ((n_haps-1)*L + trans_penalty)*np.log(N*L) - 2*model_ll
    return model, final_lls, model_ll
    
def generate_diploid_profiles(model):
    n_haps = len(model.haps)
    result = SVModel([], [])
    
    for i in range(n_haps):
        for j in range(i, n_haps):
            hap = model.haps[i] + model.haps[j]
            freq = model.freqs[i]*model.freqs[j]*(1+(i!=j))
            
            result.haps.append(hap)
            result.freqs.append(freq)
    result.normalise()
    return result

def get_sv_boundaries(fstart, fend, start, svlen, svtype, bin_size = 1000):
    if svtype == 'DEL':
        start = start
        end = start + svlen
    else:
        end = start + svlen
        start = start - svlen
        
    start_rounded = (start // bin_size) * bin_size
    end_rounded = ((end + bin_size - 1) // bin_size) * bin_size
    start_idx = (start_rounded - fstart) // bin_size
    end_idx = (end_rounded - fstart) // bin_size
    return start_idx, end_idx

def evaluate_per_hap(hap, pre_computed_lls):
    L, N = pre_computed_lls[0].shape
    
    result = np.zeros(N)
    for i in range(L):
        cn = int(hap[i])
        block = pre_computed_lls[cn]
        result = result + block[i,:]
    return result

def normalise_ll(lls, min_prob=1e-6):
    probs = np.exp(lls - logsumexp(lls, axis=1, keepdims=True))
    
    probs[probs < min_prob] = 0.0
    probs /= probs.sum(axis=1, keepdims=True)
    return probs

def get_best_haps(model, probs):
    N, K = probs.shape
    n_haps = len(model.haps)
    keys = {}
    count = 0
    for i in range(n_haps):
        for j in range(i, n_haps):
            keys[count] = (i, j)
            count += 1
    
    indices = probs.argmax(axis = 1)
    results = [keys[i] for i in indices]
    return results

def sort_model(model):
    if len(model.haps) == 1:
        return model
    else:
        haps = model.haps[1:]
        freqs = model.freqs[1:]

        sorted_haps = [hap for hap, _ in sorted(zip(haps, freqs), key=lambda x: x[1], reverse=True)]
        sorted_freqs = [freq for _, freq in sorted(zip(haps, freqs), key=lambda x: x[1], reverse=True)]
        
        new_haps = [model.haps[0]] + sorted_haps
        new_freqs = [model.freqs[0]] + sorted_freqs
    
    return SVModel(new_haps, new_freqs)

def call_sv_samples(samples, genotypes):
    results = {}
    results[(0,0)] = []
    
    for i, g in enumerate(genotypes):
        if g in results.keys():
            results[g].append(samples[i])
        else:
            results[g] = [samples[i]]
    return dict(sorted(results.items()))

def nonahore(means, variances, covs,
             ploidy = 2, 
             n_iter = 10, 
             n_sample_freq = 200, 
             n_recomb = 1000,
             bin_size = 1000,
             geom_penalty = 0.9,
             max_cnv = 10,
             min_threshold = 0.01,
             verbose = False):
    
    pre_computed_lls = precompute_site_lls(means/ploidy, variances/ploidy, covs)

    L, N = covs.shape

    reference_model = SVModel([np.ones(L)], [1])
    best_model = copy.deepcopy(reference_model)

    best_model_ary = []
    best_ll_ary = []
    best_penalty_ary = []

    for iteration in range(n_iter):
        this_sample_freq = min(10*iteration, n_sample_freq)
        this_recomb = max(n_recomb - 50*iteration, int(n_recomb/2))
        
        reference_model = copy.deepcopy(best_model)
        old_haps_set = set()
        models = [reference_model]

        if len(reference_model.haps) > 1:
            for j in range(1, len(reference_model.haps)):
                model = copy.deepcopy(reference_model)
                model.haps = model.haps[:j] + model.haps[(j+1):]
                model.freqs = model.freqs[:j] + model.freqs[(j+1):]
                model.normalise()
                models.append(model)

        if len(reference_model.haps) > 1:
            for _ in range(this_sample_freq):
                model = copy.deepcopy(reference_model)
                model.freqs = dirichlet_sampling(model.freqs)
                model.normalise()
                if sum(np.array(model.freqs) >= min_threshold) > 1:
                    models.append(model)

        for _ in range(this_recomb):
            new_hap = sample_recombinants(reference_model, L)
            if new_hap.tobytes() in old_haps_set:
                pass
            else:
                old_haps_set.add(new_hap.tobytes())
                for f in np.arange(1,20)*0.05:
                    model = copy.deepcopy(reference_model)
                    model.add(new_hap, f)
                    models.append(model)

                if len(reference_model.haps) > 1:
                    for ix in range(1, len(reference_model.haps)):
                        model = copy.deepcopy(reference_model)
                        model.replace(new_hap,ix)
                        models.append(model)

        _, _, ref_model_ll = multi_evaluate_model(reference_model, N, pre_computed_lls, geom_penalty)

        best_model = copy.deepcopy(reference_model)
        best_model_ll = ref_model_ll
        
        if sys.platform.startswith("linux"):
            ncores = max(2*(len(os.sched_getaffinity(0))) - 1, 1)
        else:
            ncores = max(2*(multiprocessing.cpu_count()-1) - 1, 1)

        with multiprocessing.Pool(processes=ncores) as pool:
            results = pool.starmap(
                multi_evaluate_model,
                [(m, N, pre_computed_lls, geom_penalty) for m in models]
            )

        for res in results:
            m, _, model_ll = res
            if model_ll < best_model_ll:
                best_model = sort_model(copy.deepcopy(m))
                best_model_ll = model_ll
                reference_model = copy.deepcopy(m)

        best_model_ary.append(reference_model)
        best_ll_ary.append(best_model_ll)
 
        if verbose and ((iteration + 1) % 5 == 0):
            print(f'------ Iteration {iteration + 1} ------')
            print(f'Best loglikelihood: {best_ll_ary[iteration]}')
            
        if len(best_ll_ary) >= 100 and len(set(best_ll_ary[-100:])) == 1:
            break
    
    _, lls, _ = multi_evaluate_model(best_model, N, pre_computed_lls, geom_penalty)
    best_probs = normalise_ll(lls)
    best_genotypes = get_best_haps(best_model, best_probs)
    
    results = {}
    results['model_ary'] = best_model_ary
    results['ll_ary'] = best_ll_ary
    results['probs'] = best_probs
    results['genotypes'] = best_genotypes
    return results

def get_ticks(df, tick_step):
    df = df.reset_index(drop = True)
    ticks = df['position'].to_numpy()
    if tick_step not in [0.001, 0.01, 0.1, 1]:
        raise ValueError(f"Invalid tick_step: {tick_step}. Choose from 0.001, 0.01, 0.1 or 1.")

    tick_start =ticks[0]
    tick_end = ticks[-1]
    ticks = ticks[np.isclose((ticks / tick_step) % 1, 0)]
    if len(ticks) == 0:
        ticks = [tick_start, tick_end]
    return ticks

def check_low_mq(arr, n = 3, score = 20):
    mask = arr < score
    window_sum = np.convolve(mask, np.ones(n, dtype=int), 'valid')
    has_consecutive = np.any(window_sum == n)
    return has_consecutive

def evaluate_sim_model(result_dict, true_hap, true_gt):
    freq = 0
    concordance = 0
    info = 0
    
    probs = result_dict['probs']
    genotypes = result_dict['genotypes']
    best_model = result_dict['model_ary'][-1]
    N = len(genotypes)
    
    if (len(best_model.haps) != 2) or all(best_model.haps[1] != true_hap):
        pass

    else:
        freq = best_model.freqs[1]
        
        dosage_ary = []
        dosage2_ary = []

        for i in range(N):
            h1, h2 = true_gt[i]
            t1, t2 = genotypes[i]
            concordance += (max(((h1 == t1) + (h2 == t2)), ((h1 == t2) + (h2 == t1))))
            
            p = probs[i]
            dosage_ary.append(p[1] + 2*p[2])
            dosage2_ary.append(p[1] + 4*p[2])
        
        concordance = concordance/N
        dosage_ary = np.array(dosage_ary)
        dosage2_ary = np.array(dosage2_ary)
        maf = np.mean(dosage_ary)/2
        info = 1 - np.mean((dosage2_ary - dosage_ary**2)/(2*maf*(1-maf)))
    return concordance, info, freq

def evaluate_real_model(result_dict):
    freq = 0
    info = 0
    
    probs = result_dict['probs']
    genotypes = result_dict['genotypes']
    best_model = result_dict['model_ary'][-1]
    N = len(genotypes)
    
    if (len(best_model.haps) != 2):
        pass
    else:
        freq = best_model.freqs[1]
        
        dosage_ary = []
        dosage2_ary = []

        for i in range(N):
            t1, t2 = genotypes[i]
            
            p = probs[i]
            dosage_ary.append(p[1] + 2*p[2])
            dosage2_ary.append(p[1] + 4*p[2])
        
        dosage_ary = np.array(dosage_ary)
        dosage2_ary = np.array(dosage2_ary)
        maf = np.mean(dosage_ary)/2
        info = 1 - np.mean((dosage2_ary - dosage_ary**2)/(2*maf*(1-maf)))
    return info, freq

def evaluate_real_model2(result_dict, plausible_boundaries, svtype):
    freqs = []
    infos = [1]
    concordance = 0
    
    probs = result_dict['probs']
    genotypes = result_dict['genotypes']
    best_model = result_dict['model_ary'][-1]
    N = len(genotypes)
    
    true_hap = None
    true_hap_idx = -9
        
    if len(best_model.haps) == 1:
        freqs = [1]
    else:
        for i, h in enumerate(best_model.haps):
            start_idx, end_idx = plausible_boundaries
            true_hap_ind = compare_plausible_sv(h, start_idx, end_idx)
            
            freqs.append(best_model.freqs[i])
            if i != 0:
                infos.append(compute_multiallelic_info(probs, i))
            
            if true_hap_ind and ((svtype == 'INS' and np.any(h >= 2)) or (svtype == 'DEL' and np.any(h == 0))):
                true_hap = h
                true_hap_idx = i
                concordance = 1
    
    return infos, freqs, concordance, true_hap_idx

def evaluate_real_model3(result_dict, plausible_boundaries, svtype, include_bins, n_bins):
    freqs = []
    infos = [1]
    concordance = 0
    
    probs = result_dict['probs']
    genotypes = result_dict['genotypes']
    best_model = result_dict['model_ary'][-1]
    N = len(genotypes)
    
    true_hap = None
    true_hap_idx = -9
        
    if len(best_model.haps) == 1:
        freqs = [1]
    else:
        for i, h in enumerate(best_model.haps):
            start_idx, end_idx = plausible_boundaries
            true_hap_ind = compare_plausible_sv2(h, start_idx, end_idx, include_bins, n_bins)
            
            freqs.append(best_model.freqs[i])
            if i != 0:
                infos.append(compute_multiallelic_info(probs, i))
            
            if check_sv2(true_hap_ind, h, svtype, plausible_boundaries, include_bins):
                true_hap = h
                true_hap_idx = i
                concordance = 1
    
    return infos, freqs, concordance, true_hap_idx

def check_sv2(true_hap_ind, h, svtype, plausible_boundaries, include_bins):
    start_idx, end_idx = plausible_boundaries
    called_bins = np.where((include_bins >= start_idx) & (include_bins < end_idx))[0].size
    c1 = true_hap_ind
    c2 = (np.where(h != 1)[0].size <= called_bins + 1)
    if svtype == 'INS':
        c2p = ((h[h > 1] - 1).sum() <= called_bins + 1)
    else:
        c2p = False
    c3 = ((svtype == 'INS' and np.any(h >= 2)) or (svtype == 'DEL' and np.any(h == 0)))
    return c1 and (c2 or c2p) and c3

def check_sv(true_hap_ind, h, svtype, plausible_boundaries, include_bins):
    start_idx, end_idx = plausible_boundaries
    called_bins = np.where((include_bins >= start_idx) & (include_bins < end_idx))[0].size
    c1 = true_hap_ind
    c2 = (np.where(h != 1)[0].size >= called_bins - 1) and (np.where(h != 1)[0].size <= called_bins + 1)
    if svtype == 'INS':
        c2p = ((h[h > 1] - 1).sum() >= called_bins - 1) and ((h[h > 1] - 1).sum() <= called_bins + 1)
    else:
        c2p = False
    c3 = ((svtype == 'INS' and np.any(h >= 2)) or (svtype == 'DEL' and np.any(h == 0)))
    
    return c1 and (c2 or c2p) and c3

def compare_plausible_sv2(arr, start_idx, end_idx, include_bins, n_bins):
    full = np.full(n_bins, None)
    for i, b in enumerate(include_bins):
        full[b] = arr[i]

    sub_bins = np.arange(start_idx, end_idx)
    sub_vals = [full[b] for b in sub_bins if full[b] is not None]

    if not sub_vals:
        return False

    non1 = np.array(sub_vals) != 1
    padded = np.r_[False, non1, False]
    starts = np.where(~padded[:-1] & padded[1:])[0]
    ends   = np.where(padded[:-1] & ~padded[1:])[0]

    return len(starts) == 1
#     return (len(starts) == 1) and (np.unique(arr).size == 2)

def compare_plausible_sv(arr, start_idx, end_idx):
    sub = arr[start_idx:end_idx]
    non1 = np.array(sub) != 1
    padded = np.r_[False, non1, False]
    starts = np.where(~padded[:-1] & padded[1:])[0]
    ends = np.where(padded[:-1] & ~padded[1:])[0]
    return len(starts) == 1
#     return (len(starts) == 1) and (np.unique(arr).size == 2)

def compute_multiallelic_info(GP, alt_index):
    N, K = GP.shape
    A = int((np.sqrt(8 * K + 1) - 1) / 2)
    if A * (A + 1) // 2 != K:
        raise ValueError("Invalid number of genotype columns in GP.")
    if alt_index < 1 or alt_index >= A:
        raise ValueError(f"alt_index must be in [1, {A-1}] for {A} alleles.")

    genotypes = list(combinations_with_replacement(range(A), 2))
    alt_copy_per_genotype = np.array([g.count(alt_index) for g in genotypes])

    dosage = GP @ alt_copy_per_genotype
    dosage2 = GP @ (alt_copy_per_genotype ** 2)

    var_g = np.mean(dosage2 - dosage**2)
    p = np.mean(dosage) / 2

    if p == 0 or p == 1:
        return 0.0

    info = 1 - var_g / (2 * p * (1 - p))
    return info

def estimate_maf(means, variances, coverage, svtype, fillna = True):
    this_vars = coverage.var(axis = 1)

    seq_means = means/2
    seq_vars = variances/2

    EM = seq_means.mean()
    DM = seq_means.var(ddof = 1)
    V = seq_vars.mean()

    A = 2*DM - 2*(EM**2)
    
    if svtype == 'DEL':
        B = 2*(EM**2-V-4*DM+DM)
    else:
        B = 2*(EM**2+V+4*DM+DM)
    C = 4*DM + 2*V-this_vars

    delta = B**2-4*A*C
    sqrt_term = np.sqrt(delta)
    f1 = (-B + sqrt_term) / (2*A)
    f2 = (-B - sqrt_term) / (2*A)
    
    if fillna:
        f1 = np.nan_to_num(f1, nan=0.0)
    
    return f1

def find_intervals(means, f, min_run=5):
    crit1 = (means >= means.mean() / 2)
    crit2 = (f >= 0.05)
    crit = (crit1 & crit2).astype(np.uint8)

    arr = np.asarray(crit, dtype=bool)
    padded = np.concatenate(([False], arr, [False]))

    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    lengths = ends - starts
    
    valid = lengths >= min_run
    valid_starts = starts[valid]
    valid_lengths = lengths[valid]

    df = pd.DataFrame({
        'POS': valid_starts,
        'SVLEN': valid_lengths
    })   
    return df

def find_overlapping_svs(df, c, start, end):
    gr = pr.PyRanges(df)
    gr_chr = gr[gr.Chromosome == c]
    match = (start, end)

    input_row = df[(df["Chromosome"] == c) &
                    (df["Start"] == start) &
                    (df["End"] == end)]

    df_clustered = gr_chr.cluster().df
    query_cluster_id = df_clustered[
        (df_clustered["Start"] == start) &
        (df_clustered["End"] == end)
    ]["Cluster"].values[0]

    cluster_df = df_clustered[df_clustered["Cluster"] == query_cluster_id]
    return cluster_df["Start"].min(), cluster_df["End"].max()

def print_avg_concordance(df):
    print('DEL:', df[df['SVTYPE'] == 'DEL']['concordance'].mean())
    print('INS:', df[df['SVTYPE'] == 'INS']['concordance'].mean())
    return None

def plot_sv_coverage(cov, chromosome, start, end, flank, calling_dict, side = 'both', tick_step = 0.1):
    means, _ = normalise_by_flank(cov, start, end, flank, side = side)
    
    region = cov[(cov['position'] >= start) & (cov['position'] <= end)]
    region['position'] = region['position']/1e6
    
    colors = plt.get_cmap(CATEGORY_CMAP_STR).colors[:10]
    colors = [mcolors.to_hex(c) for c in colors]
    
    all_samples = np.array([s.split(':')[0] for s in cov.columns[1:-1]])

    for i, k in enumerate(calling_dict.keys()):
        samples = calling_dict[k]
        for s in samples:
            index = np.where(all_samples == s)[0][0]
            
            if k == (0,0):
                plt.plot(region['position'], region[f'{s}:coverage']/means[index], alpha = 1, color = '0.8')
            else:
                plt.plot(region['position'], region[f'{s}:coverage']/means[index], alpha = 1, color = colors[i - 1])
   
    ticks = get_ticks(region, tick_step)
    if tick_step == 0.001:
        plt.xticks(ticks, [f"{int(tick*1000)}" for tick in ticks], rotation = 45)
        plt.xlabel(f'Chromosome {chromosome} (Kb)')
    else:
        plt.xticks(ticks, [f"{tick:.{int(np.log10(1/tick_step))}f}" for tick in ticks], rotation = 45)
        plt.xlabel(f'Chromosome {chromosome} (Mb)')

    color_handles = []
    color_index = [0,0]
    for i, k in enumerate(calling_dict.keys()):
        if k != (0,0):
            color_handles.append(Line2D([0], [0], color=colors[i-1], label=k))

    legend1 = plt.legend(handles=color_handles, loc='upper left', prop={'size': 10}, framealpha=1)
    legend1.get_title().set_fontsize(9)
    plt.gca().add_artist(legend1)

    plt.ylabel('Coverage (X)')
    return None

def plot_sv_coverage_by_gt(cov, chromosome, start, end, flank, calling_dict, side = 'both', tick_step = 0.1):
    means, _ = normalise_by_flank(cov, start, end, flank, side = side)
    
    region = cov[(cov['position'] >= start) & (cov['position'] <= end)]
    region['position'] = region['position']/1e6
    region.iloc[:,1:-1] = region.iloc[:,1:-1]/means[np.newaxis,:]
    
    if len(calling_dict.keys()) > 10:
        print('Only region with less than 3 different haplotypes can be printed.')
        return None
    
    colors = plt.get_cmap(CATEGORY_CMAP_STR).colors[:10]
    colors = [mcolors.to_hex(c) for c in colors]
    
    all_samples = np.array([s.split(':')[0] for s in region.columns[1:-1]])
    for i, k in enumerate(calling_dict.keys()):
        samples = calling_dict[k]
        indices = np.where(np.isin(all_samples, samples))[0] + 1
        
        tmp = pd.concat([region.iloc[:, 0], region.iloc[:,indices]], axis = 1)
        plt.plot(tmp['position'], tmp.iloc[:,1:].mean(axis = 1), alpha = 1, color = colors[i])
   
    ticks = get_ticks(region, tick_step)
    if tick_step == 0.001:
        plt.xticks(ticks, [f"{int(tick*1000)}" for tick in ticks], rotation = 45)
        plt.xlabel(f'Chromosome {chromosome} (Kb)')
    else:
        plt.xticks(ticks, [f"{tick:.{int(np.log10(1/tick_step))}f}" for tick in ticks], rotation = 45)
        plt.xlabel(f'Chromosome {chromosome} (Mb)')

    color_handles = []
    color_index = [0,0]
    for i, k in enumerate(calling_dict.keys()):
        color_handles.append(Line2D([0], [0], color=colors[i], label=k))

    legend1 = plt.legend(handles=color_handles, loc='upper left', prop={'size': 10}, framealpha=1)
    legend1.get_title().set_fontsize(9)
    plt.gca().add_artist(legend1)

    plt.ylabel('Coverage (X)')
    return None

def plot_training(results, show_legends = True):
    L = len(results['model_ary'][0].haps[0])
    lls = results['ll_ary']
    n = len(lls)
    haps = {np.ones(L).tobytes(): 0}
    freqs = np.zeros((1, n))

    for i in range(n):
        m = results['model_ary'][i]
        for j, hap in enumerate(m.haps):
            h = hap.tobytes()
            if h not in haps.keys():
                haps[h] = freqs.shape[0]
                freqs = np.append(freqs, np.zeros((1,n)), axis=0)
            ridx = haps[h]
            freqs[ridx, i] = m.freqs[j]

    n_haps = freqs.shape[0]
    x = np.arange(1, n+1)
    fig, ax1 = plt.subplots()

    ax1.plot(x, lls, ls = '--', color='black')
    ax1.set_xlabel('Iterations')
    ax1.set_ylabel('Model score')
    ax1.tick_params(axis='y')
    y1_max = lls[0]
    y1_min = lls[-1]
    y1_ext = (y1_max - y1_min)/3
    ax1.set_ylim((y1_min - y1_ext, y1_max + y1_ext))

    colors = plt.get_cmap(CATEGORY_CMAP_STR).colors[:n_haps]
    colors = [mcolors.to_hex(c) for c in colors]

    ax2 = ax1.twinx()
    for i in range(freqs.shape[0]):
        x = np.flatnonzero(freqs[i,:])
        ax2.plot(x, freqs[i,x]*100, color=colors[i])
    ax2.set_ylabel('Haplotype frequencies (%)')
    ax2.tick_params(axis='y')
    ax2.set_ylim((-5, 105))

    ax1.grid(True, alpha = 0.7)
    
    if show_legends:
        color_handles = []
        for i in range(n_haps):
            color_handles.append(Line2D([0], [0], color=colors[i], label=f'Haplotype {i}'))

        legend1 = plt.legend(handles=color_handles, title='Haplotypes', 
                             prop={'size': 10}, framealpha=1)
        legend1.get_title().set_fontsize(10)
        plt.gca().add_artist(legend1)

        linestyle_handles = [
            Line2D([0], [0], color='black', lw=2, linestyle='--', label='Model score')
        ]
        legend2 = plt.legend(handles=linestyle_handles, title='Modified BIC', 
                             prop={'size': 10}, framealpha=1, bbox_to_anchor = (1, 0.8))
        legend2.get_title().set_fontsize(10)
        plt.gca().add_artist(legend2)   
        legend2.get_frame().set_zorder(2)
    return None

def plot_sv_heatmap(means, variances, coverage, samples, results, n_sample = 10):
    best_model = results['model_ary'][-1]
    L, N = coverage.shape
    pre_computed_lls = precompute_site_lls(means, variances, coverage)
    calling_dict = call_sv_samples(samples, results['genotypes'])
    s1 = calling_dict[(0,0)][:n_sample]
    s2 = calling_dict[(0,1)][:n_sample]

    best_model = results['model_ary'][-1]
    n_haps = len(best_model.haps)
    L = len(best_model.haps[0])
    n_diploid = int(n_haps*(n_haps + 1)/2)

    dip_model = generate_diploid_profiles(best_model)

    g1 = np.zeros((int(n_sample*2), L))
    for i, s in enumerate(s1 + s2):
        all_samples = np.array(samples)
        sample_index = np.where(all_samples == s)[0][0]

        h1 = dip_model.haps[0]
        for l in range(L):
            cn = int(h1[l])
            block = pre_computed_lls[cn]
            g1[i, l] = block[l, sample_index]

    ax = sns.heatmap(
        g1,
        cmap='magma',
        xticklabels=False,
        yticklabels=True,
        cbar_kws={'label': 'Log-likelihood'}
    )
    ax.axhline(y=n_sample, color='black', linewidth=2)
    plt.xlabel("Chromosomal bins")
    plt.ylabel("Individuals")
    plt.tight_layout()
    
    return None

def plot_maf(coverage, f1):
    this_vars = coverage.var(axis = 1) 
    fig, ax1 = plt.subplots()

    ax1.plot(this_vars, color = 'b', label = 'VarC')
    ax1.legend()
    ax1.set_ylabel('Realised variance of coverage')
    ax2 = ax1.twinx()
    ax2.plot(f1, color = 'orange', label = 'MAF est.')
    ax2.legend()
    ax2.set_ylabel('Estimated SV allele frequency')
    
    return None

def plot_sv_result_category(tmp, info_cutoff = 0.8):
    counts = pd.DataFrame(np.zeros((2,4)),columns = ['Monomorphic', 'Wrong SV', f'True SV (INFO<{info_cutoff})', f'True SV (INFO>{info_cutoff})'])

    Ns = []
    for i, t in enumerate(['DEL', 'INS']):
        tmp1 = tmp[tmp['SVTYPE'] == t]
        N = len(tmp1)
        Ns.append(N)

        counts.iloc[i, 0] = (tmp1['n_hap'] == 1).sum()/N
        counts.iloc[i, 1] = (tmp1['concordance'] == 0).sum()/N - counts.iloc[i, 0]

        tmp2 = tmp1[tmp1['concordance'] != 0]
        counts.iloc[i, 2] = (tmp2['true_info'] <= info_cutoff).sum()/N
        counts.iloc[i, 3] = (tmp2['true_info'] > info_cutoff).sum()/N

    rgba_colors = plt.cm.RdYlBu(np.linspace(0, 1, 4))
    hex_codes = [mcolors.to_hex(c) for c in rgba_colors]

    ax = counts.plot(kind='bar', stacked=True, color = hex_codes, alpha = 0.7)
    for bar_idx, row in enumerate(counts.values):
        cumulative = 0
        for col_idx, value in enumerate(row):
            percent = value*100
            if value > 0:  
                ax.text(
                    bar_idx,                    
                    cumulative + value / 2,
                    f"{percent:.1f}%",
                    ha='center', va='center', color='black', fontsize=10
                )
            cumulative += value

    plt.xticks(np.arange(2), [f'Deletion\nN={Ns[0]}', f'Insertion\nN={Ns[1]}'], rotation = 0)
    plt.yticks(np.arange(6)*0.2, np.arange(6)*20)
    plt.ylabel('Proportion (%)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    return None

def plot_mappability(df):
    x = []
    y = []
    for _, row in df.iterrows():
        x.extend([row["start"], row["end"]])
        y.extend([row["value"], row["value"]])

    plt.figure(figsize=(8,3))
    plt.plot(x, y, drawstyle="steps-post")
    plt.ylim(-0.05, 1.05)
    plt.xlabel("Genomic position")
    plt.ylabel("Mappability")
    return None


def calculate_score_per_alignment(cigars, cg = -6, cm = 0, cx = -4, ce = -2):
    score = 0
    
    for i, c in enumerate(cigars):
        code, l = c
        if code == 1 or code == 2 or code == 4:
            score += ce*l
            score += cg
        elif code == 8:
            score += cx*l
        else:
            pass
    return score  