import io
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import csv
from scipy.stats import poisson
import itertools
import collections

def get_imputed_dosage(df, colname='call'):
    ref = df['ref']
    alt = df['alt']
    s = df[colname]
    if alt == '.' or len(alt) > 1 or len(ref) > 1:
        return pd.NaT  
    if s[:3] == '0|0' or s[:3] == '0/0':
        return pd.NaT
    else:
        return s.split(':')[2]
    
def get_genotype(df, colname='call'):
    ref = df['ref']
    alt = df['alt']
    s = df[colname]
    if alt == '.' or len(alt) > 1 or len(ref) > 1:
        return pd.NaT
    if s[:3] == '0|0' or s[:3] == '0/0':
        return pd.NaT
    elif s[:3] == '1|0' or s[:3] == '1/0':
        return 1
    elif s[:3] == '0|1' or s[:3] == '0/1':
        return 1
    elif s[:3] == '1|1' or s[:3] == '1/1':
        return 2
    else:
        return pd.NaT
    
def get_imputed_dosage_72175193(df):
    return get_imputed_dosage(df, '72175193')
def get_imputed_dosage_72245200(df):
    return get_imputed_dosage(df, '72245200')
def get_imputed_dosage_72605236(df):
    return get_imputed_dosage(df, '72605236')
def get_imputed_dosage_72655241(df):
    return get_imputed_dosage(df, '72655241')
def get_genotype_72135285(df):
    return get_genotype(df, '72135285')
def get_genotype_72125284(df):
    return get_genotype(df, '72125284')
def get_genotype_72115283(df):
    return get_genotype(df, '72115283')
def get_genotype_72105282(df):
    return get_genotype(df, '72105282')

for i in range(1,23):
    tmp = pd.read_csv('results/variant_calling/lcwgs_1x_calling/lcwgs_1x_variant_calling_chr' + str(i) + '.txt', header = None, sep = '\t', names = ['chr', 'pos', 'ref', 'alt', 'attributes', '72175193', '72245200', '72605236', '72655241'],
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
        })
    tmp['72175193'] = tmp.apply(get_imputed_dosage_72175193, axis=1)
    tmp['72245200'] = tmp.apply(get_imputed_dosage_72245200, axis=1)
    tmp['72605236'] = tmp.apply(get_imputed_dosage_72605236, axis=1)
    tmp['72655241'] = tmp.apply(get_imputed_dosage_72655241, axis=1)
    tmp = tmp.dropna()
    tmp['72175193'] = pd.to_numeric(tmp['72175193'])
    tmp['72245200'] = pd.to_numeric(tmp['72245200'])
    tmp['72605236'] = pd.to_numeric(tmp['72605236'])
    tmp['72655241'] = pd.to_numeric(tmp['72655241'])
    tmp = tmp.drop(columns = ['attributes'])
    if i==1:
        lcwgs_1x_vcf = tmp
    else:
        lcwgs_1x_vcf = pd.concat([lcwgs_1x_vcf,tmp])

for i in range(1,23):
    tmp = pd.read_csv('results/variant_calling/lcwgs_20x_calling/lcwgs_20x_snps_chr'+str(i)+'.txt', header = None, sep = '\t', 
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
    tmp['72105282'] = pd.to_numeric(tmp['72105282'])
    tmp['72115283'] = pd.to_numeric(tmp['72115283'])
    tmp['72125284'] = pd.to_numeric(tmp['72125284'])
    tmp['72135285'] = pd.to_numeric(tmp['72135285'])
    if i==1:
        lcwgs_20x_vcf = tmp
    else:
        lcwgs_20x_vcf = pd.concat([lcwgs_20x_vcf,tmp])

for i in range(1,23):
    tmp = pd.read_csv('results/variant_calling/gnomAD_MAFs/gnomAD_MAF_chr'+str(i)+'.txt', header = None, sep = '\t', names = ['chr', 'pos', 'ref', 'alt', 'MAF'],
                      dtype = {
        'chr': 'string',
        'pos': 'Int64',
        'ref': 'string',
        'alt': 'string',
        'MAF': 'string'
    })
    tmp = tmp.dropna()
    tmp['MAF'] = pd.to_numeric(tmp['MAF'])
    tmp = tmp.dropna()
    if i==1:
        gnomAD_vcf = tmp
    else:
        gnomAD_vcf = pd.concat([gnomAD_vcf,tmp])

lcwgs = pd.merge(lcwgs_1x_vcf, lcwgs_20x_vcf, on=['chr', 'pos', 'ref', 'alt'], how="left")
lcwgs = pd.merge(lcwgs, gnomAD_vcf, on=['chr', 'pos', 'ref', 'alt'], how="left")
lcwgs = lcwgs.dropna()

MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])
r2 = np.zeros((4, np.size(MAF_ary) - 1))
for i in range(r2.shape[1]):
    tmp = lcwgs[(MAF_ary[i+1] > lcwgs['MAF']) & (lcwgs['MAF'] > MAF_ary[i])]
    r2[0,i] = np.corrcoef(tmp['72105282'].values, tmp['72175193'].values)[0,1]**2
    r2[1,i] = np.corrcoef(tmp['72115283'].values, tmp['72245200'].values)[0,1]**2
    r2[2,i] = np.corrcoef(tmp['72125284'].values, tmp['72605236'].values)[0,1]**2
    r2[3,i] = np.corrcoef(tmp['72135285'].values, tmp['72655241'].values)[0,1]**2

plt.plot(MAF_ary[1:], r2[0,:], label = '72175193', color = 'c')
plt.plot(MAF_ary[1:], r2[1,:], label = '72245200', color = 'g')
plt.plot(MAF_ary[1:], r2[2,:], label = '72605236', color = 'r')
plt.plot(MAF_ary[1:], r2[3,:], label = '72655241', color = 'm')
plt.xlabel('gnomAD MAF')
plt.ylabel('$r^2$')
plt.legend()
plt.title('Imputation Accuracy of 4 samples')
plt.xscale('log')
plt.savefig('graphs/lcwgs_imputation_accuracy.png', bbox_inches = "tight", dpi=300)