import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import csv
from scipy.stats import poisson
import itertools
import collections

def get_genotype(df, colname='call'):
    ref = df['ref']
    alt = df['alt']
    s = df[colname]
    if alt == '.':
        alt = ref
    if s[:3] == '0|0' or s[:3] == '0/0':
        return ref+'_'+ref
    elif s[:3] == '1|0' or s[:3] == '1/0':
        return ref+'_'+alt
    elif s[:3] == '0|1' or s[:3] == '0/1':
        return ref+'_'+alt
    elif s[:3] == '1|1' or s[:3] == '1/1':
        return alt+'_'+alt
    else:
        return pd.NaT
def get_genotype_72175193(df):
    return get_genotype(df, '72175193')
def get_genotype_72245200(df):
    return get_genotype(df, '72245200')
def get_genotype_72605236(df):
    return get_genotype(df, '72605236')
def get_genotype_72655241(df):
    return get_genotype(df, '72655241')
def get_genotype_72135285(df):
    return get_genotype(df, '72135285')
def get_genotype_72125284(df):
    return get_genotype(df, '72125284')
def get_genotype_72115283(df):
    return get_genotype(df, '72115283')
def get_genotype_72105282(df):
    return get_genotype(df, '72105282')
def get_genotype_10x(df):
    return get_genotype(df, '10x')
def get_genotype_15x(df):
    return get_genotype(df, '15x')

# lcwgs_1x
i = 0
for chunk in pd.read_csv('variant_calling/lcwgs_1x_calling/lcwgs_1x.txt', header = None, sep = '\t', 
                            names = ['chr', 'pos', 'ref', 'alt', 'attributes', '72175193', '72245200', '72605236', '72655241'],
                      dtype = {
            'chr': 'string',
            'pos': 'Int64',
            'ref': 'string',
            'alt': 'string',
            'attributes': 'string',
            '72175193': 'string',
            '72245200': 'string',
            '72605236': 'string',
            '72655241': 'string'
        }, chunksize=1e6):
    chunk['EAF'] = chunk.attributes.str.extract('EAF=([^;]+)')
    chunk['score'] = chunk.attributes.str.extract('INFO_SCORE=([^;]+)')
    chunk = chunk.drop(columns = ['attributes'])
    chunk['72175193_genotype'] = chunk.apply(get_genotype_72175193, axis=1)
    chunk['72245200_genotype'] = chunk.apply(get_genotype_72245200, axis=1)
    chunk['72605236_genotype'] = chunk.apply(get_genotype_72605236, axis=1)
    chunk['72655241_genotype'] = chunk.apply(get_genotype_72655241, axis=1)
    chunk = chunk.dropna()
    chunk = chunk.drop(columns = ['72175193', '72245200','72605236','72655241', 'ref', 'alt'])
    chunk = chunk.iloc[:,[0,1,4,5,6,7,2,3]]
    chunk.iloc[:,6] = pd.to_numeric(chunk.iloc[:,6])
    chunk.iloc[:,7] = pd.to_numeric(chunk.iloc[:,7])
    chunk.rename(columns={'72605236_genotype': '72605236', '72655241_genotype': '72655241', '72175193_genotype': '72175193', '72245200_genotype': '72245200'}, inplace=True)
    if i == 0:
        lcwgs_1x_vcf = chunk
        i += 1
    else:
        lcwgs_1x_vcf = pd.concat([lcwgs_1x_vcf, chunk])

# lcwgs_20x
for i in range(1,23):
    tmp = pd.read_csv('variant_calling/lcwgs_20x_calling/lcwgs_20x_snps_chr'+str(i)+'.txt', header = None, sep = '\t', 
                      names = ['chr', 'pos', 'ref', 'alt', '72105282', '72115283', '72125284', '72135285'],
                      dtype = {
            'chr': 'string',
            'pos': 'Int64',
            'ref': 'string',
            'alt': 'string',
            '72105282': 'string',
            '72115283': 'string',
            '72125284': 'string',
            '72135285': 'string'
        })
    tmp['72105282'] = tmp.apply(get_genotype_72105282, axis=1)
    tmp['72115283'] = tmp.apply(get_genotype_72115283, axis=1)
    tmp['72125284'] = tmp.apply(get_genotype_72125284, axis=1)
    tmp['72135285'] = tmp.apply(get_genotype_72135285, axis=1)
    tmp = tmp.dropna()
    tmp = tmp.drop(columns = ['ref', 'alt'])
    if i==1:
        lcwgs_20x_vcf = tmp
    else:
        lcwgs_20x_vcf = pd.concat([lcwgs_20x_vcf,tmp])

lcwgs_vcf = pd.merge(lcwgs_1x_vcf, lcwgs_20x_vcf, on=['chr', 'pos'], how="left", suffixes=('_1x', '_20x'))
lcwgs_vcf = lcwgs_vcf.dropna()

def check_maf(maf):
    if maf<0.0005:
        return 0
    elif maf<0.001:
        return 1
    elif maf<0.005:
        return 2
    elif maf<0.01:
        return 3
    elif maf<0.05:
        return 4
    elif maf<0.1:
        return 5
    elif maf<0.5:
        return 6
    else:
        return 7

def calculate_imputation_accuracy(df):
    match_mat = np.zeros((4,8))
    total_mat = np.zeros((4,8))
    for i in range(df.shape[0]):
        ix = check_maf(df.iloc[i,6])
        total_mat[:, ix] += 1
        for j in range(4):
            if df.iloc[i,2+j] == df.iloc[i,8+j]:
                 match_mat[j, ix] += 1
    prop_mat = match_mat/total_mat
    return prop_mat

prop_mat = calculate_imputation_accuracy(lcwgs_vcf)

fig = plt.figure(figsize=(6,4))
EAF_ary = np.array([0, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5])
plt.plot(EAF_ary, prop_mat[0,:]*100, label = '72175193', color = 'c')
plt.plot(EAF_ary, prop_mat[1,:]*100, label = '72245200', color = 'g')
plt.plot(EAF_ary, prop_mat[2,:]*100, label = '72605236', color = 'r')
plt.plot(EAF_ary, prop_mat[3,:]*100, label = '72655241', color = 'm')
plt.xlabel('MAF')
plt.ylabel('Imputation Accuracy (%)')
plt.legend()
plt.title('Imputation Accuracy')
plt.xscale('log')
plt.figtext(.5, -0.05, 
            'Figure 8. Imputation accuracy of lcwgs samples.', ha='center', wrap = True)
plt.savefig('imputation_accuracy_lcwgs.png', bbox_inches = "tight")
