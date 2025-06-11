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

def read_coverage_data(file_path, sep = ','):    
    df = pd.read_csv(file_path, sep = sep)
    bins = list(zip(df['position'], df['position'] + df['size'], df['N']))
    samples = [col for col in df.columns if col.endswith(":coverage") and col != "total:coverage"]
    coverage = df.loc[:,samples].replace("NA", 0).astype(float).to_numpy()
    samples = [s.split(':')[0] for s in samples]
    return samples, bins, coverage

def extract_target_cov(df, start, end):
    df = df[(df['position'] >= start) & (df['position'] <= end)].reset_index(drop = True)
    samples = [col for col in df.columns if col.endswith(":coverage") and col != "total:coverage"]
    coverage = df.loc[:,samples].replace("NA", 0).astype(float).to_numpy()
    samples = [s.split(':')[0] for s in samples]
    return samples, coverage

def normalise_by_flank(df, start, end, flank, side = 'both'):
    fstart = max(start-flank, df.iloc[0,0])
    fend = min(end+flank, df.iloc[-1,0])
    
    if side == 'both':
        criteria = (((df['position'] >= fstart) & (df['position'] < start)) | 
             ((df['position'] > end) & (df['position'] <= fend)))
    elif side == 'left':
        criteria = ((df['position'] >= fstart) & (df['position'] < start))
    elif side == 'right':
        criteria = ((df['position'] > end) & (df['position'] <= fend))
    else:
        raise ValueError("Unsupported side.")            
            
    cov = df[criteria].iloc[:,1:-1].to_numpy()
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

def precompute_site_lls(means, variances, coverage, max_cnv=10, mismap_proportion = 0.01, model='negative_binomial'):
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
        
def dirichlet_sampling(freqs, concentration = 10):
    alpha = np.array(freqs) * concentration
    rng = default_rng()
    return list(rng.dirichlet(alpha))

def sample_from_freqs(freqs):
    freqs = np.array(freqs)
    cum_probs = np.cumsum(freqs)
    u = np.random.uniform(0, 1)
    return np.searchsorted(cum_probs, u)

def sample_recombinants(model, L, max_cnv = 10):
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
    return new_hap

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
    return results

def nonahore(means, variances, covs,
             ploidy = 2, 
             n_iter = 10, 
             n_sample_freq = 200, 
             n_recomb = 1000,
             bin_size = 1000,
             geom_penalty = 0.9,
             max_cnv = 10,
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
            
        if len(best_ll_ary) >= 50 and len(set(best_ll_ary[-50:])) == 1:
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
    plt.xticks(ticks, [f"{tick:.{int(np.log10(1/tick_step))}f}" for tick in ticks], rotation = 45)

    color_handles = []
    color_index = [0,0]
    for i, k in enumerate(calling_dict.keys()):
        if k != (0,0):
            color_handles.append(Line2D([0], [0], color=colors[i-1], label=k))

    legend1 = plt.legend(handles=color_handles, loc='upper left', prop={'size': 10}, framealpha=1)
    legend1.get_title().set_fontsize(9)
    plt.gca().add_artist(legend1)

    plt.xlabel(f'Chromosome {chromosome} (Mb)')
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
    ax1.set_ylabel('logL')
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

        legend1 = plt.legend(handles=color_handles, title='Haplotype frequencies', 
                             prop={'size': 10}, framealpha=1)
        legend1.get_title().set_fontsize(10)
        plt.gca().add_artist(legend1)

        linestyle_handles = [
            Line2D([0], [0], color='black', lw=2, linestyle='--', label='logL')
        ]
        legend2 = plt.legend(handles=linestyle_handles, title='Training Process', 
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