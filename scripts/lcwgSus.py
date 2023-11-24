import io
import os
import sys
import csv
import gzip
import time
import random
import secrets
import resource
import itertools
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import statsmodels.api as sm
import scipy
from scipy.stats import poisson
from scipy.stats import chi2
from scipy.stats import friedmanchisquare
from scipy.stats import studentized_range
pd.options.mode.chained_assignment = None

def get_mem():
    current_memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    current_memory_usage_mb = current_memory_usage / 1024
    print(f"Current memory usage: {current_memory_usage_mb:.2f} MB")
def get_genotype(df, colname='call'):
    ref = df['ref']
    alt = df['alt']
    s = df[colname]
    if len(alt) != 1 or len(ref) != 1:
        return np.NaN
    if s == '0|0' or s == '0/0':
        return 0.
    elif s == '1|0' or s == '1/0':
        return 1.
    elif s == '0|1' or s == '0/1':
        return 1.
    elif s == '1|1' or s == '1/1':
        return 2.
    else:
        return np.NaN
def get_imputed_dosage(df, colname='call'): # Currently deprecated
    ref = df['ref']
    alt = df['alt']
    s = df[colname]
    if alt == '.' or len(alt) > 1 or len(ref) > 1 :
        return pd.NaT
    else:
        return s.split(':')[2]
def convert_to_str(x):
    if x == int(x):
        return str(int(x))
    else:
        return str(x)
def read_metadata(file, filetype = 'gzip', comment = '#'):
    if filetype == 'gzip':
        with io.TextIOWrapper(gzip.open(file,'r')) as f:
            metadata = [l for l in f if l.startswith(comment)]
    else:
        with open(file, 'r') as f:
            metadata = [l for l in f if l.startswith(comment)]
    return metadata
def read_vcf(file, sample = 'call', q = None): 
    colname = read_metadata(file)
    header = colname[-1].replace('\n', '').split('\t')
    df = pd.read_csv(file, compression='gzip', comment='#', sep = '\t', header = None, names = header).rename(columns={'#CHROM': 'chr', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt'})
#     with io.TextIOWrapper(gzip.open(file,'r')) as f:
#         lines =[l for l in f if not l.startswith('##')]
#         dynamic_header_as_key = []
#         for liness in f:
#             if liness.startswith("#CHROM"):
#                 dynamic_header_as_key.append(liness)
#         values = [str,int,str,str,str,int,str,str,str,str]
#         columns2detype = dict(zip(dynamic_header_as_key,values))
#         df = pd.read_csv(
#             io.StringIO(''.join(lines)),
#             dtype=columns2detype,
#             sep='\t'
#         ).rename(columns={'#CHROM':'CHROM'})
    if df.dtypes[0] != int: # Continue for now, but need to change this later if we are not merely considering autosomes
        df['CHROM'] = df['CHROM'].str.extract(r'(\d+)').astype(int)
    if df.dtypes[1] != int:
        df['POS'] = df['POS'].astype(int)
    if len(df.columns) == 10: 
        df.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'call']
        if sample != 'call':
            df.columns[-1] = sample
    if q is None:
        return df
    else:
        q.put(df)
def save_vcf(df, metadata, save_name = 'test.vcf.gz'): # Only use this if no cols are removed from the original vcf
    # df is the vcf_df to be saved
    # metadata is a list generated from read_metadata
    if type(df.iloc[0,0] != str):
        df[df.columns[0]] = 'chr' + df[df.columns[0]].astype(str)
    random_str = secrets.token_hex(8) + '_'
    file_path = random_str + 'test.vcf'
    metadata_path = random_str + 'metadata.txt'
    df.to_csv(file_path, index=False, sep = '\t', header = False)
    with open(metadata_path, 'w') as metadata_file:
        metadata_file.write(''.join(metadata))
    gzipped_file_path = save_name
    with open(metadata_path, 'rb') as metadata_file:
        metadata_content = metadata_file.read()
    with open(file_path, 'rb') as data_file, gzip.open(gzipped_file_path, 'wb') as gzipped_file:
        gzipped_file.write(metadata_content)
        gzipped_file.writelines(data_file)
    os.remove(file_path)
    os.remove(metadata_path)
def read_af(file, q = None):
    df = pd.read_csv(file, header = None, sep = '\t', names = ['chr', 'pos', 'ref', 'alt', 'MAF'],
                      dtype = {
        'chr': 'string',
        'pos': 'Int64',
        'ref': 'string',
        'alt': 'string',
        'MAF': 'string'
    })
    df = df.dropna()
    df['MAF'] = pd.to_numeric(df['MAF'])
    df['chr'] = df['chr'].str.extract(r'(\d+)').astype(int)
    if q is None:
        return df
    else:
        q.put(df)
def read_r2(panels, samples, drop=3): # Modify indir
    dfs = []
    for i in panels:
        for j in samples:
            tmp = pd.read_csv("results/imputation/imputation_accuracy/"+j+"/"+i+"_imputation_accuracy.csv", sep = ',', dtype = {
                'MAF': float,
                'Imputation Accuracy': float,
                'Bin Count': str
            }).iloc[drop:,:]
            tmp['panel'] = i
            tmp['Bin Count'] = j
            tmp.columns = ['AF', 'corr', 'sample', 'panel']
            tmp['AF'] = (100*tmp['AF']).apply(convert_to_str)
            tmp['AF'] = tmp['AF'].shift(1).fillna('0') + '-' + tmp['AF']
            tmp['AF'] = tmp['AF'].astype("category")
            dfs.append(tmp)
    res = pd.concat(dfs).reset_index(drop = True)
    return res
def aggregate_r2(df):
    tmp = df.copy().groupby(['AF', 'panel'])['corr'].mean().reset_index()
    res_ary = []
    for i in tmp['panel'].unique():
        imp_res = tmp[tmp['panel'] == i]
        imp_res['sort'] = imp_res['AF'].apply(lambda x: x.split('-')[0]).astype(float)
        imp_res = imp_res.sort_values(by = 'sort', ascending = True).drop(columns = 'sort')
        res_ary.append(imp_res)
    return res_ary
def extract_info(df, info_cols = ['EAF', 'INFO_SCORE'], attribute = 'info', drop_attribute = True):
    for i in info_cols:
        df[i] = df[attribute].str.extract( i + '=([^;]+)' ).astype(float)
    if drop_attribute:
        df = df.drop(columns = [attribute])
    return df
def extract_format(df, sample, fmt = 'format'):
    fields = df[fmt].values[0].split(':')
    try:
        df[fields] = df[sample].str.split(':', expand=True)
        df[df.columns[-1]] = df[df.columns[-1]].astype(float)
        if len(fields) != len(df[sample].values[0].split(':')):
            raise ValueError("Mismatching fields in FORMAT and Imputed results.")
    except ValueError as e:
        print(f"Error: {e}")
    return df.drop(columns = [fmt, sample])
def drop_cols(df, drop_lst = ['id', 'qual', 'filter']):
    return df.drop(columns = drop_lst)

def parse_vcf(file, sample = 'call', q = None, 
              info_cols = ['EAF', 'INFO_SCORE'], attribute = 'info', fmt = 'format', drop_attribute = True, drop_lst = ['id', 'qual', 'filter']):
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

def multi_parse_vcf(chromosomes, files, parse = True, sample = 'call', combine = True,
               info_cols = ['EAF', 'INFO_SCORE'], attribute = 'info', fmt = 'format', drop_attribute = True, drop_lst = ['id', 'qual', 'filter']):
    manager = multiprocessing.Manager()
    q = manager.Queue()
    processes = []
    for i in range(len(chromosomes)):
        if parse:
            tmp = multiprocessing.Process(target=parse_vcf, args=(files[i], sample, q, info_cols, attribute, fmt, drop_attribute, drop_lst))
        else:
            tmp = multiprocessing.Process(target=read_vcf, args=(files[i], sample, q))
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
def multi_read_af(chromosomes, files, combine = True):
    manager = multiprocessing.Manager()
    q = manager.Queue()
    processes = []
    for i in range(len(chromosomes)):
        tmp = multiprocessing.Process(target=read_af, args=(files[i], q))
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
def calculate_af(df, drop = True):
    # df should have columns chr, pos, ref, alt and genotypes
    df['prop'] = 0
    for i in range(len(df.index)):
        count = 0
        for j in range(4, len(df.columns) - 2):
            count += df.iloc[i, j].split('/').count('1')
        df.iloc[i, -1] = count/(2*(len(df.columns) - 5))
    if drop:
        return df[['chr', 'pos', 'ref', 'alt', 'prop']]
    else:
        return df
def filter_afs(df1, df2, diff = 0.2, z_score = None):
    # df1 is the main vcf in which afs are to be filtered out
    # df2 is the ref panel afs
    # Either filter by z-score (suggested 2 sds so 1.96 or diff=0.2)
    res = pd.merge(df1, df2, on = ['chr', 'pos', 'ref', 'alt'])
    if z_score is not None:
        res = res[(res['prop_y'] != 0) & (res['prop_y'] != 1)]
        res['z'] = (res['prop_x'] - res['prop_y'])/np.sqrt(res['prop_y']*(1-res['prop_y']))
        res = res[abs(res['z']) <= z_score]
        return res.drop(columns = ['prop_x', 'prop_y', 'z'])
    else:
        res = res[abs(res['prop_x'] - res['prop_y']) < diff]
        return res.drop(columns = ['prop_y']).rename(columns = {'prop_x': 'prop'})
def calculate_ss_cumsum_coverage(df, num_coverage=5):
    df['bases'] = df['end'] - df['start']
    df = df.groupby(['cov']).bases.sum().reset_index()
    df['prop bases'] = df['bases']/df.bases.sum()
    df['cum prop'] = np.cumsum(df['prop bases'].to_numpy())
    df['prop genome at least covered'] = (1-df['cum prop'].shift(1))
    df = df.dropna()
    coverage_ary = df['prop genome at least covered'].values[:num_coverage]
    return coverage_ary
def calculate_imputation_accuracy(df1, df2, af,
                                  MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.95, 1]),
                                 how = 'left'):
    df2 = df2.copy()
    if len(df1.columns) != 5:
        df1 = df1[['chr', 'pos', 'ref', 'alt', 'DS']]
    col1 = df1.columns[-1]
    if type(df2.iloc[0, len(df2.columns)-1]) == str:
        df2['genotype'] = df2.apply(get_genotype, axis = 1)
        df2 = df2.dropna()
        df2['genotype'] = df2['genotype'].astype(float)
        df2 = df2.drop(columns = df2.columns[-2])
        col2 = 'genotype'
    else:
        col2 = df2.columns[-1]

    df = pd.merge(df2, df1, on=['chr', 'pos', 'ref', 'alt'], how=how)
    df = df.fillna(0)
    df = pd.merge(df, af, on=['chr', 'pos', 'ref', 'alt'], how='left')
    df = df.dropna()

    r2 = np.zeros((2, np.size(MAF_ary) - 1))
    for i in range(r2.shape[1]):
        tmp = df[(MAF_ary[i+1] > df['MAF']) & (df['MAF'] > MAF_ary[i])]
        if tmp.shape[0] == 0:
            r2[0,i] = 0
        else:
            r2[0, i] = np.corrcoef(tmp[col1].values, tmp[col2].values)[0,1]**2
        r2[1, i] = int(tmp.shape[0])

    r2_df = pd.DataFrame(r2.T, columns = ['Imputation Accuracy','Bin Count'], index = MAF_ary[1:])
    r2_df.index.name = 'MAF'
    return r2_df
def plot_afs(df1, df2, save_fig = False, outdir = 'graphs/', save_name = 'af_vs_af.png'):
    # df1 is the chip df with cols chr, pos, ref, alt and prop
    # df2 is the other df with the same cols
    df = pd.merge(df1, df2, on = ['chr', 'pos', 'ref', 'alt'], how = 'inner')
    plt.scatter(df['prop_x']*100, df['prop_y']*100)
    plt.xlabel('ChIP MAF (%)')
    plt.ylabel('GGVP AF (%)')
    plt.title('Check AFs')
    if save_fig:
        plt.savefig(outdir + save_name, bbox_inches = "tight", dpi=300)
    return np.corrcoef(df['prop_x'], df['prop_y'])[0,1]
def plot_imputation_accuracy(r2, single_sample = True, aggregate = True, save_fig = False, save_name = 'imputation_corr_vs_af.png', outdir = 'graphs/'):
    plt.figure(figsize = (10,6))
    if single_sample:
        if type(r2) == pd.DataFrame:
            plt.plot(r2.index, r2['Imputation Accuracy'], color = 'g')
        else:
            for i in range(len(r2)):
                plt.plot(r2[i].index, r2[i]['Imputation Accuracy'])
        plt.xlabel('gnomAD AF (%)')
        plt.ylabel('$r^2$')
        plt.title(plot_title)
        plt.xscale('log')
    else:
        if aggregate:
            for df in r2:
                panel = df['panel'].values[0]
                plt.plot(np.arange(1, df.shape[0]+1), df['corr'], label = panel)
            plt.xticks(np.arange(1, r2[0].shape[0]+1), r2[0]['AF'], rotation = 45)
            plt.xlabel('Allele frequencies (%)')
            plt.legend()
            plt.text(x = -1.5, y = 1.02, s = 'Aggregated imputation accuracy ($r^2$)')
            plt.grid(alpha = 0.5)
        else:
            sns.set(style="whitegrid")
            sns.stripplot(data=r2, x="corr", y="AF", hue="panel", dodge=True)
            plt.xlabel('Imputation Accuracy')
            plt.ylabel('gnomAD allele frequencies')
    if save_fig:
        plt.savefig(outdir + save_name, bbox_inches = "tight", dpi=300)
def plot_sequencing_skew(arys, avg_coverage, n_se = 1.96, code = None, num_coverage=5, save_fig = False, save_name = 'prop_genome_at_least_coverage.png', outdir = 'graphs/'):
    poisson_expectation = 1 - np.cumsum(poisson.pmf(np.arange(num_coverage), mu=avg_coverage, loc=0))
    se = np.sqrt(avg_coverage/len(arys))
    x_coordinate = np.arange(1, num_coverage+1)
    plt.figure(figsize=(16,12))
    for i in range(len(arys)):
        coverage_ary = arys[i]
        plt.plot(x_coordinate, coverage_ary/poisson_expectation[0], label = code) # Can put code in as well
    plt.plot(x_coordinate, poisson_expectation/poisson_expectation[0], label = 'Poisson', ls='--', color = 'k', linewidth = 5)
    plt.plot(x_coordinate, (poisson_expectation + n_se*se)/poisson_expectation[0], ls='--', color = 'k', linewidth = 5)
    plt.plot(x_coordinate, (poisson_expectation - n_se*se)/poisson_expectation[0], ls='--', color = 'k', linewidth = 5)
    plt.xticks(x_coordinate)
    plt.xlabel('Coverage (x)')
    plt.ylabel('Sequencing Skew')
    #plt.legend()
    plt.title('Sequencing Skew')
    if save_fig:
        plt.savefig(outdir + save_name, bbox_inches = "tight", dpi=300)
def plot_info_vs_af(vcf, afs, MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.95, 1]),
                   save_fig = False, outdir = 'graphs/', save_name = 'info_vs_af.png'):
    df = pd.merge(vcf[['chr', 'pos', 'ref', 'alt', 'info']], afs, on=['chr', 'pos', 'ref', 'alt'], how="left").dropna()
    df['classes'] = np.digitize(df['MAF'], MAF_ary)
    plt.figure(figsize = (12,8))
    sns.violinplot(data=df, x="classes", y="info")
    plt.xlabel('Allele Frequencies (%)')
    plt.ylabel('INFO_SCORE')
    plt.title('INFO Score vs Allele Frequencies')
    ax = plt.gca()
    ax.set_xticklabels(MAF_ary[np.sort(df['classes'].unique()) - 1])
    if save_fig:
        plt.savefig(outdir + save_name, bbox_inches = "tight", dpi=300)
