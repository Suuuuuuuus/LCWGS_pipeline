import io
import os
import csv
import gzip
import time
import multiprocessing
import resource
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import random
from collections import Counter
import matplotlib.colors as mcolors
from scipy.stats import poisson
import itertools
import collections
import scipy
from scipy.stats import chi2
from scipy.stats import friedmanchisquare
from scipy.stats import studentized_range
pd.options.mode.chained_assignment = None

def get_mem():
    current_memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    current_memory_usage_mb = current_memory_usage / 1024
    print(f"Current memory usage: {current_memory_usage_mb:.2f} MB")
def read_vcf(file):
    with io.TextIOWrapper(gzip.open(file,'r')) as f:
        lines =[l for l in f if not l.startswith('##')]
        dynamic_header_as_key = []
        for liness in f:
            if liness.startswith("#CHROM"):
                dynamic_header_as_key.append(liness)
        values = [str,int,str,str,str,int,str,str,str,str]
        columns2detype = dict(zip(dynamic_header_as_key,values))
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype=columns2detype,
            sep='\t'
        ).rename(columns={'#CHROM':'CHROM'})
    df['CHROM'] = df['CHROM'].str.extract(r'(\d+)').astype(int)
    return df
def extract_info(df, info_cols = ['EAF', 'INFO_SCORE'], attribute = 'INFO', drop_attribute = True):
    for i in info_cols:
        df[i] = df[attribute].str.extract( i + '=([^;]+)' ).astype(float)
    if drop_attribute:
        df = df.drop(columns = [attribute])
    return df
def extract_format(df, sample, fmt = 'FORMAT'):
    fields = df[fmt].values[0].split(':')
    try:
        df[fields] = df[sample].str.split(':', expand=True)
        df[df.columns[-1]] = df[df.columns[-1]].astype(float)
        if len(fields) != len(df[sample].values[0].split(':')):
            raise ValueError("Mismatching fields in FORMAT and Imputed results.")
    except ValueError as e:
        print(f"Error: {e}")
    return df.drop(columns = [fmt, sample])
def drop_cols(df, drop_lst = ['ID', 'QUAL', 'FILTER']):
    return df.drop(columns = drop_lst)

def parse_vcf(file, sample, q = None, 
              info_cols = ['EAF', 'INFO_SCORE'], attribute = 'INFO', fmt = 'FORMAT', drop_attribute = True, drop_lst = ['ID', 'QUAL', 'FILTER']):
    df = read_vcf(file)
    df = extract_info(df, info_cols = info_cols, attribute = attribute, drop_attribute = drop_attribute)
    df = extract_format(df, sample, fmt = fmt)
    df = drop_cols(df, drop_lst = drop_lst)
    if q is None:
        return df
    else:
        q.put(df)
def file_to_list(df):
    lst = []
    for i in df[df.columns[0]].unique():
        lst.append(df[df[df.columns[0]] == i])
    return lst
def combine_df(lst):
    df = lst[0]
    for i in range(1, len(lst)):
        df = pd.concat([df, lst[i]])
    return df.sort_values(by = df.columns[:2].to_list()).reset_index(drop = True)
def subtract_bed_by_chr(cov, region, q = None):
    i = 0
    tmp = 0
    for j in range(region.shape[0]):
        chr, start, end = region.iloc[j,:]
        while start > cov.iloc[i,2]:
            i += 1
        if start < cov.iloc[i,1]:
            cov.iloc[i-1, 2] = start
            if end < cov.iloc[i,2]:
                cov.iloc[i,1] = end
            elif end == cov.iloc[i,2]:
                cov.iloc[i,3] = -9
                i += 1
            else:
                tmp = i
                while end > cov.iloc[tmp,2]:
                    tmp += 1
                if end < cov.iloc[tmp, 2]:
                    cov.iloc[tmp, 1] = end
                    cov.iloc[i:tmp, 3] = -9
                    i = tmp
                else: 
                    cov.iloc[i:tmp+1, 3] = -9
                    i = tmp
        elif start == cov.iloc[i,1]:
            if end < cov.iloc[i,2]:
                cov.iloc[i,1] = end
            elif end == cov.iloc[i,2]:
                cov.iloc[i, 3] = -9
            else:
                tmp = i
                while end > cov.iloc[tmp,2]:
                    tmp += 1
                if end < cov.iloc[tmp, 2]:
                    cov.iloc[tmp, 1] = end
                    cov.iloc[i:tmp+1, 3] = -9
                    i = tmp
                else: 
                    cov.iloc[i:tmp, 3] = -9
                    i = tmp
        else: 
            idx = cov.index.max() + 1
            cov.loc[idx] = {'chr': chr, 'start': cov.iloc[i,1], 'end': start, 'cov': cov.iloc[i,3]}
            if end < cov.iloc[i, 2]:
                cov.iloc[i, 1] = end
            elif end == cov.iloc[i, 2]:
                cov.iloc[i, 3] = -9
            else: 
                tmp = i
                while end > cov.iloc[tmp,2]:
                    tmp += 1
                if end < cov.iloc[tmp, 2]:
                    cov.iloc[tmp, 1] = end
                    cov.iloc[i:tmp, 3] = -9
                    i = tmp
                else: 
                    cov.iloc[i:tmp+1, 3] = -9
                    i = tmp
    res = cov[cov['cov'] >= 0].sort_values(by = cov.columns[:2].to_list()).reset_index(drop = True)
    if q is None:
        return res
    else:
        q.put(res)

def multi_parse(chromosomes, files, sample, combine = True,
               info_cols = ['EAF', 'INFO_SCORE'], attribute = 'INFO', fmt = 'FORMAT', drop_attribute = True, drop_lst = ['ID', 'QUAL', 'FILTER']):
    manager = multiprocessing.Manager()
    q = manager.Queue()
    processes = []
    for i in range(len(chromosomes)):
        tmp = multiprocessing.Process(target=parse_vcf, args=(files[i], sample, q, info_cols, attribute, fmt, drop_attribute, drop_lst))
        tmp.start()
        processes.append(tmp)
    for process in processes:
        process.join()
    res_lst = []
    while not q.empty():
        res_lst.append(q.get())
    if combine:
        return combine_df(res_lst)
    else:
        return res_lst
def multi_subtract_bed(chromosomes, covs, regions, combine = True):
    manager = multiprocessing.Manager()
    q = manager.Queue()
    processes = []
    for i in range(len(chromosomes)):
        tmp = multiprocessing.Process(target=subtract_bed_by_chr, args=(covs[i], regions[i], q))
        tmp.start()
        processes.append(tmp)
    for process in processes:
        process.join()
    res_lst = []
    while not q.empty():
        res_lst.append(q.get())
    if combine:
        return combine_df(res_lst)
    else:
        return res_lst
    
def calculate_ss_cumsum_coverage(df, num_coverage=5):
    df['bases'] = df['end'] - df['start']
    df = df.groupby(['cov']).bases.sum().reset_index()
    df['prop bases'] = df['bases']/df.bases.sum()
    df['cum prop'] = np.cumsum(df['prop bases'].to_numpy())
    df['prop genome at least covered'] = (1-df['cum prop'].shift(1))
    df = df.dropna()
    coverage_ary = df['prop genome at least covered'].values[:num_coverage]
    return coverage_ary
def plot_depth_coverage(arys, avg_coverage, n_se = 1.96, code = None, num_coverage=5, save_fig = False, outdir = 'graphs/'):
    poisson_expectation = 1 - np.cumsum(poisson.pmf(np.arange(num_coverage), mu=avg_coverage, loc=0))
    se = np.sqrt(avg_coverage/len(arys))
    x_coordinate = np.arange(1, num_coverage+1)
    plt.figure(figsize=(16,12))
    for i in range(len(arys)):
        coverage_ary = arys[i]
        plt.plot(x_coordinate, coverage_ary*100, label = code) # Can put code in as well
    plt.plot(x_coordinate, poisson_expectation*100, label = 'Poisson', ls='--', color = 'k', linewidth = 5)
    plt.plot(x_coordinate, (poisson_expectation + n_se*se)*100, ls='--', color = 'k', linewidth = 5)
    plt.plot(x_coordinate, (poisson_expectation - n_se*se)*100, ls='--', color = 'k', linewidth = 5)
    plt.xticks(x_coordinate)
    plt.xlabel('Coverage (x)')
    plt.ylabel('Proportion of genome at least coverage (%)')
    #plt.legend()
    plt.title('Genome coverage')
    if save_fig:
        plt.savefig(outdir + 'prop_genome_at_least_coverage.png', bbox_inches = "tight", dpi=300)