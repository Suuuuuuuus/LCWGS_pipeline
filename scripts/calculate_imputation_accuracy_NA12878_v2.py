import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
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
        return 0
    elif s[:3] == '1|0' or s[:3] == '1/0':
        return 1
    elif s[:3] == '0|1' or s[:3] == '0/1':
        return 1
    elif s[:3] == '1|1' or s[:3] == '1/1':
        return 2
    else:
        return pd.NaT

for i in range(1,23):
    tmp = pd.read_csv('results/variant_calling/NA12878_1x_calling/72155287_variant_calling_chr' + str(i) + '.txt', header = None, sep = '\t', names = ['chr', 'pos', 'ref', 'alt', 'attributes', 'call'],
                          dtype = {
            'chr': 'string',
            'pos': 'Int64',
            'ref': 'string',
            'alt': 'string',
            'attributes': 'string',
            'call': 'string'
        })
    tmp['dosage'] = tmp.apply(get_imputed_dosage, axis=1)
    tmp = tmp.dropna()
    tmp['dosage'] = pd.to_numeric(tmp['dosage'])
    tmp = tmp.drop(columns = ['call', 'attributes'])
    if i==1:
        NA12878_1x_vcf = tmp
    else:
        NA12878_1x_vcf = pd.concat([NA12878_1x_vcf,tmp])

NA12878_20x_vcf = pd.read_csv('results/variant_calling/genome_in_the_bottle/NA12878_all.txt', header = None, sep = '\t', names = ['chr', 'pos', 'ref', 'alt', 'call'],
                    dtype = {
    'chr': 'string',
    'pos': 'Int64',
    'ref': 'string',
    'alt': 'string',
    'call': 'string'
})
NA12878_20x_vcf['call'] = NA12878_20x_vcf.apply(get_genotype, axis=1)
NA12878_20x_vcf = NA12878_20x_vcf.dropna()
NA12878_20x_vcf['call'] = pd.to_numeric(NA12878_20x_vcf['call'])

NA12878_high_confidence = pd.read_csv('results/variant_calling/genome_in_the_bottle/NA12878_high_confidence_position.txt', header = None, sep = '\t', names = ['chr', 'pos'],
                      dtype = {'chr': 'string', 'pos': 'Int64'})

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
    if i==1:
        gnomAD_vcf = tmp
    else:
        gnomAD_vcf = pd.concat([gnomAD_vcf,tmp])

NA12878 = pd.merge(NA12878_1x_vcf, NA12878_20x_vcf, on=['chr', 'pos', 'ref', 'alt'], how="left")
NA12878 = pd.merge(NA12878, gnomAD_vcf, on=['chr', 'pos', 'ref', 'alt'], how="left")
NA12878 = NA12878.dropna()
NA12878_high_confidence = pd.merge(NA12878_high_confidence, NA12878, on=['chr', 'pos'], how="left")
NA12878_high_confidence = NA12878_high_confidence.dropna()

MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])
r2 = np.zeros(np.size(MAF_ary) - 1)
r2_high_confidence = np.zeros(np.size(MAF_ary) - 1)
for i in range(r2.size):
    tmp = NA12878[(MAF_ary[i+1] > NA12878['MAF']) & (NA12878['MAF'] > MAF_ary[i])]
    r2[i] = np.corrcoef(tmp['dosage'].values, tmp['call'].values)[0,1]**2
    tmp = NA12878_high_confidence[(MAF_ary[i+1] > NA12878_high_confidence['MAF']) & (NA12878_high_confidence['MAF'] > MAF_ary[i])]
    r2_high_confidence[i] = np.corrcoef(tmp['dosage'].values, tmp['call'].values)[0,1]**2

plt.plot(MAF_ary[1:], r2, label = 'all', color = 'g')
plt.plot(MAF_ary[1:], r2_high_confidence, label = 'high_confidence', color = 'r')
plt.xlabel('gnomAD MAF')
plt.ylabel('$r^2$')
plt.legend()
plt.title('Imputation Accuracy of the NA12878 sample')
plt.xscale('log')
plt.savefig('graphs/NA12878_imputation_accuracy.png', bbox_inches = "tight", dpi = 300)
