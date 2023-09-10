import io
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import random
from collections import Counter
import csv
import glob
import matplotlib.colors as mcolors
from scipy.stats import poisson
import itertools
import collections

path = sys.argv[1]
num_coverage = int(sys.argv[2])
avg_coverage = float(sys.argv[3])

rstrip = '_subsampled_bedgraph.txt'
lstrip = 'results/coverage/subsampled_bedgraphs/'
code = path.replace(rstrip, "").replace(lstrip, "")

df = pd.read_csv(path , header = None, sep = '\t', 
                 names = ['chr', 'start', 'end', 'cov'],
                dtype = {
                    'chr': 'string',
                    'start': 'Int64',
                    'end': 'Int64',
                    'cov': 'Int64'
                })
df = df[df.chr.isin([str(i) for i in range(1,23)] + ['chr' + str(i) for i in range(1,23)])]
df['bases'] = df['end'] - df['start']
df = df.groupby(['cov']).bases.sum().reset_index()
df['prop bases'] = df['bases']/df.bases.sum()
df['cum prop'] = np.cumsum(df['prop bases'].to_numpy())
df['prop genome at least covered'] = (1-df['cum prop'].shift(1))
df = df.dropna()
coverage_ary = df['prop genome at least covered'].values[:num_coverage]
np.savetxt('results/coverage/subsampled_bedgraphs/' + code + '_cumsum_ary.txt', coverage_ary, newline = ' ', fmt='%1.7f')
