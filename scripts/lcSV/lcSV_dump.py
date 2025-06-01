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