import io
import os
import sys
import csv
import gzip
import time
import json
import secrets
import copy
import multiprocessing
import subprocess
import resource
import itertools
import collections
import sqlite3
import random

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

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


best_model = results['model_ary'][-1]
L, N = coverage.shape
pre_computed_lls = precompute_site_lls(means, variances, coverage)
s1 = calling_dict[(0,0)][:10]
s2 = calling_dict[(0,1)][:10]

best_model = results['model_ary'][-1]
n_haps = len(best_model.haps)
L = len(best_model.haps[0])
n_diploid = int(n_haps*(n_haps + 1)/2)

dip_model = generate_diploid_profiles(best_model)

g1 = np.zeros((20, L))
g2 = np.zeros((20, L))
for i, s in enumerate(s1):
    all_samples = np.array(samples)
    sample_index = np.where(all_samples == s)[0][0]
    
    h1 = dip_model.haps[0]
    h2 = dip_model.haps[1]
    for l in range(L):
        cn = int(h1[l])
        block = pre_computed_lls[cn]
        g1[i, l] = block[l, sample_index]
        
        cn = int(h2[l])
        block = pre_computed_lls[cn]
        g2[i, l] = block[l, sample_index]

for i, s in enumerate(s2):
    all_samples = np.array(samples)
    sample_index = np.where(all_samples == s)[0][0]
    
    h1 = dip_model.haps[1]
    h2 = dip_model.haps[0]
    for l in range(L):
        cn = int(h1[l])
        block = pre_computed_lls[cn]
        g1[i+10, l] = block[l, sample_index]
        
        cn = int(h2[l])
        block = pre_computed_lls[cn]
        g2[i+10, l] = block[l, sample_index]

sns.heatmap(
    g1,
    cmap='YlGnBu',
    xticklabels=True,  # Hide tick labels if too many bins
    yticklabels=True,
    cbar_kws={'label': 'Log-likelihood'}
)

plt.xlabel("Genomic Bin")
plt.ylabel("Individual")
plt.title("Log-Likelihood Heatmap per Individual Across Genomic Bins")
plt.tight_layout()
plt.show()

sns.heatmap(
    g2,
    cmap='YlGnBu',
    xticklabels=True,  # Hide tick labels if too many bins
    yticklabels=True,
    cbar_kws={'label': 'Log-likelihood'}
)

plt.xlabel("Genomic Bin")
plt.ylabel("Individual")
plt.title("Log-Likelihood Heatmap per Individual Across Genomic Bins")
plt.tight_layout()
plt.show()


def evaluate_model(model, N, pre_computed_lls, geom_penalty, bin_size = 1000):
    n_haps = len(model.haps)
    L = len(model.haps[0])
    n_diploid = int(n_haps*(n_haps + 1)/2)
    
    dip_model = generate_diploid_profiles(model)
    final_lls = np.zeros((N, n_diploid))
    for i in range(n_diploid):
        ll_ary = evaluate_per_hap(dip_model.haps[i], pre_computed_lls)
        final_lls[:, i] = ll_ary

    final_lls = final_lls + np.log(np.array(dip_model.freqs))[np.newaxis, :] 
    model_ll = logsumexp(final_lls, axis=1).sum()
#     model_ll = model_ll + np.log(geom.pmf(n_haps, geom_penalty))
    model_ll = ((n_haps-1)*L + trans)*np.log(N*L) - 2*model_ll
    return final_lls, model_ll

def main(regions, sv_df_file, eichler_file, ix, ofile, bin_size = 1000):
    eichler_full = read_eichler(eichler_file)
    indices = np.where(eichler_full['SVTYPE'] == 'INS')[0]
    eichler_full.loc[indices, 'END'] = eichler_full.loc[indices, 'POS'] + eichler_full.loc[indices, 'SVLEN']
    eichler_full.loc[indices, 'POS'] = eichler_full.loc[indices, 'POS'] - eichler_full.loc[indices, 'SVLEN']
    cols = ['Chromosome', 'Start', 'End']
    eichler_full.columns = cols + eichler_full.columns[3:].tolist()
    eichler_full = eichler_full.sort_values(by = cols).reset_index(drop = True)
    
    manifest = pd.read_csv(sv_df_file, sep = '\t')
    chrom, sv_start, L, svtype = manifest.loc[ix, ['#CHROM', 'POS', 'SVLEN', 'SVTYPE']]
    chromosome = int(chrom.replace('chr', ''))
    start, end = delineate_region2(eichler_full, chrom, sv_start, L, svtype)
    flank = 1e6
    plausible_boundaries = get_sv_boundaries(start, end, sv_start, L, svtype)

    cov = load_region_files(regions, chromosome, start, end, flank = flank)
    
    tmp = cov[(cov['position'] >= start) & (cov['position'] <= end)]
    mq = tmp['total:mean_mq'].to_numpy()
    mean_mq = tmp['total:mean_mq'].mean()
    
    cov = cov[['position'] + list(cov.columns[cov.columns.str.contains('coverage')])]
    means, variances = normalise_by_flank2(cov, start, end, flank)
    samples, coverage = extract_target_cov(cov, start, end)

    results = nonahore(means, variances, coverage, n_recomb = 1000, n_iter = 2000, verbose = False)

    probs, genotypes = results['probs'], results['genotypes']
    info, freq, concordance, true_hap_idx = evaluate_real_model2(results, plausible_boundaries, svtype)
    haps = results['model_ary'][-1].haps

    outputs = {}
    outputs['chromosome'] = chrom
    outputs['g_start'] = sv_start
    outputs['start'] = start
    outputs['end'] = end
    outputs['mq'] = mq
    outputs['mean_mq'] = mean_mq
    outputs['length'] = L
    outputs['svtype'] = svtype
    outputs['means'] = means
    outputs['variances'] = variances
    outputs['coverage'] = coverage
    outputs['info'] = info
    outputs['freq'] = freq
    outputs['concordance'] = concordance
    outputs['true_hap_idx'] = true_hap_idx
    outputs['haps'] = haps
    outputs['probs'] = np.round(probs, 4)
    outputs['genotypes'] = genotypes

    with open(ofile, 'wb') as of:
        pickle.dump(outputs, of)
    
if __name__ == "__main__":
    chunk_file = snakemake.params.chunk_file
    with open(chunk_file, 'r') as file:
        regions = json.load(file)

    odir = snakemake.params.odir
    os.makedirs(odir,  exist_ok=True)

    row_ix = int(snakemake.params.row_ix)
    sv_df_file = snakemake.input.sv_df_file
    eichler_file = snakemake.params.eichler_file
    
    ofile = snakemake.output.pickle

    main(regions, sv_df_file, eichler_file, row_ix, ofile)