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

num_coverage = int(sys.argv[1])
avg_coverage = float(sys.argv[2])

poisson_expectation = 1 - np.cumsum(poisson.pmf(np.arange(0, num_coverage), mu=avg_coverage, loc=0))

def calculate_diff_squared(coverage, poisson=poisson_expectation):
    if coverage.size != poisson.size:
        poisson = poisson[:coverage.size]
    return np.sum((coverage-poisson)**2)

cov_dict = {}

path = "results/coverage/subsampled_bedgraphs/*_bedgraph.txt"
for filename in glob.glob(path):
    rstrip = '_subsampled_bedgraph.txt'
    lstrip = 'results/coverage/subsampled_bedgraphs/'
    code = filename.replace(rstrip, "").replace(lstrip, "")

    df = pd.read_csv(filename , header = None, sep = '\t', 
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
    r2 = calculate_diff_squared(coverage_ary)
    cov_dict[code] = [coverage_ary, r2]

values_list = list(cov_dict.values())
keys_list = list(cov_dict.keys())
last_elems = [lst[-1] for lst in cov_dict.values()]
max_r2 = max(last_elems)
min_r2 = min(last_elems)
max_r2_index = last_elems.index(max_r2)
min_r2_index = last_elems.index(min_r2)
max_code = keys_list[max_r2_index]
min_code = keys_list[min_r2_index]
max_ary = values_list[max_r2_index][0]
min_ary = values_list[min_r2_index][0]

x_coordinate = np.arange(1,num_coverage+1)

plt.figure(figsize=(6,4))
plt.plot(np.arange(1,max_ary.size+1), max_ary*100, label = max_code)
plt.plot(np.arange(1,min_ary.size+1), min_ary*100, label = min_code)
plt.plot(x_coordinate, poisson_expectation*100, label = 'Poisson', ls='--')
plt.xticks(x_coordinate)
plt.xlabel('Coverage (x)')
plt.ylabel('Proportion of genome at least coverage (%)')
plt.legend()
plt.title('Genome coverage')
plt.figtext(.5, -0.1, 
            'Figure 8. Proportion of the genome that is at least x-axis covered.\nLargest r2: '+str(round(max_r2, 4))+' ; Smallest r2: '+str(round(min_r2, 4))+'.', ha='center', wrap = True)
plt.savefig('graphs/fig8_prop_genome_at_least_coverage.png', bbox_inches = "tight", dpi=300)
