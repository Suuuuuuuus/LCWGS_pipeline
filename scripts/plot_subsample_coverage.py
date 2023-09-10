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

cov_dict = {}

path = "results/coverage/subsampled_bedgraphs/*_cumsum_ary.txt"
for filename in glob.glob(path):
    rstrip = '_cumsum_ary.txt'
    lstrip = 'results/coverage/subsampled_bedgraphs/'
    code = filename.replace(rstrip, "").replace(lstrip, "")
    cov_dict[code] = np.loadtxt(filename)

plt.figure(figsize=(16,12))
for code, coverage_ary in cov_dict.items():
    x_coordinate = np.arange(1, coverage_ary.size+1)
    plt.plot(x_coordinate, coverage_ary*100, label = code)
x_coordinate = np.arange(1, num_coverage+1)
plt.plot(x_coordinate, poisson_expectation*100, label = 'Poisson', ls='--')
plt.xticks(x_coordinate)
plt.xlabel('Coverage (x)')
plt.ylabel('Proportion of genome at least coverage (%)')
#plt.legend()
plt.title('Genome coverage')
plt.figtext(.5, -0.1,
            'Figure 8. Proportion of the genome that is at least x-axis covered.', ha='center', wrap = True)
plt.savefig('graphs/fig8_prop_genome_at_least_coverage.png', bbox_inches = "tight", dpi=300)
