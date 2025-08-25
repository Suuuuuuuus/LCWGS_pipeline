import io
import os
import sys
import csv
import gzip
import time
import json
import secrets
import pickle
import multiprocessing
import subprocess
import resource
import itertools
from itertools import combinations_with_replacement
import pyranges as pr
import collections
import sqlite3
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import scipy as sp
import pandas as pd
import statsmodels.api as sm
import random
from collections import Counter
from collections import defaultdict
import copy
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

from scipy.stats import nbinom, norm, geom, beta, poisson
from scipy.special import logsumexp

sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/lcwgsus/')
import lcwgsus
from lcwgsus.variables import *

sys.path.append('/well/band/users/rbx225/GAMCC/scripts/lcSV/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/Low Coverage Data/gamcc/scripts/lcSV/')
from lcSV import *

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

pd.options.mode.chained_assignment = None

def main(regions, sv_df_file, eichler_file, ix, ofile, bin_size = 1000):
    eichler_full = parse_eichler_full()
    
    manifest = pd.read_csv(sv_df_file, sep = '\t')
    chrom, sv_start, L, svtype = manifest.loc[ix, ['#CHROM', 'POS', 'SVLEN', 'SVTYPE']]
    chromosome = int(chrom.replace('chr', ''))
    start, end = delineate_region2(eichler_full, chrom, sv_start, L, svtype)
    flank = 1e6
    plausible_boundaries = get_sv_boundaries(start, end, sv_start, L, svtype)
    start_idx, end_idx = plausible_boundaries

    mapp = '/well/band/users/rbx225/recyclable_files/eichler_sv/mappability.bed.gz'

    tsv = read_bed_tabix(mapp, chromosome, start, end)
    metrics = calculate_mappability_per_bin(tsv, start, end)

    cov = load_region_files(regions, chromosome, start, end, flank = flank)
    metrics = pd.merge(metrics, cov[['position', 'total:mean_mq']], on = ['position']).fillna(0)
    
    tmp = metrics.loc[start_idx:end_idx,:]
    callable_windows = np.where((tmp['mappability'] >= 0.9) & (tmp['total:mean_mq'] >= 20))[0].size/len(tmp)
    
    to_call, include_bins = assess_bins(metrics, plausible_boundaries)
    n_bins = len(metrics)

    tmp = cov[(cov['position'] >= start) & (cov['position'] <= end)]
    mq = tmp['total:mean_mq'].to_numpy()
    mean_mq = mq[start_idx:end_idx].mean()
    
    cov = cov[['position'] + list(cov.columns[cov.columns.str.contains('coverage')])]
    flanking_regions = cov[(cov['position'] >= start - flank) & (cov['position'] <= end + flank)]
    flanking_regions['mean'] = flanking_regions.iloc[:,1:-1].mean(axis = 1)
    flanking_regions = flanking_regions[['position', 'mean']]
    f1 = flanking_regions[(flanking_regions['position'] >= start) & (flanking_regions['position'] < end)]['mean'].mean()
    f2 = flanking_regions[(flanking_regions['position'] < start) | (flanking_regions['position'] >= end)]['mean'].mean()
    mean_ratio = f1/f2

    outputs = {}
    outputs['chromosome'] = chrom
    outputs['g_start'] = sv_start
    outputs['start'] = start
    outputs['end'] = end
    outputs['mq'] = mq
    outputs['mean_mq'] = mean_mq
    outputs['length'] = L
    outputs['svtype'] = svtype
    outputs['to_call'] = to_call
    outputs['callable_windows'] = callable_windows
    outputs['mean_ratio'] = mean_ratio

    # Remove this line below!
    to_call = True
    if (not to_call) or (include_bins.size <= 4):
        print('Nonahore will not call this region..')
        with open(ofile, 'wb') as of:
            pickle.dump(outputs, of)
        return
    
    means, variances = normalise_by_flank2(cov, start, end, flank)
    samples, coverage = extract_target_cov(cov, start, end)
    coverage = coverage[include_bins, :]
    results = nonahore(means, variances, coverage, n_recomb = 1000, n_iter = 2000, verbose = False)

    probs, genotypes = results['probs'], results['genotypes']
    info, freq, concordance, true_hap_idx = evaluate_real_model3(results, plausible_boundaries, svtype, include_bins, n_bins)
    haps = results['model_ary'][-1].haps

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
