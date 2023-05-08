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

chr = sys.argv[1]
repeat_mask_bin_size = int(sys.argv[2])
ylim = float(sys.argv[3])

cov_dict = {}

path = "results/coverage/per_bin_coverage/1x/*chr" + chr + "_base.txt"
for filename in glob.glob(path):
    rstrip = '_chr' + chr +'_base.txt'
    code = filename.replace(rstrip, "").replace('results/coverage/per_bin_coverage/1x\\', '')
    cov_dict[code] = np.loadtxt(filename)

repeat_mask = pd.read_csv('data/repeat_mask_region_by_chromosome/repeat_mask_chr' + chr + '_region.txt', header = None, names = ['start', 'end'], sep = '\t',
                         dtype = {
    'start': 'Int64',
    'end': 'Int64'
})
repeat_mask['length'] = repeat_mask['end'] - repeat_mask['start']
repeat_mask = repeat_mask[repeat_mask['length'] >= repeat_mask_bin_size]

coordinate_ary = np.loadtxt(glob.glob("results/coverage/per_bin_coverage/1x/*chr" + chr + "_coordinate.txt")[0])

fig, ax = plt.subplots(2, 1, figsize = (12,8), gridspec_kw={'height_ratios': [5, 1]})
for code, coverage_ary in cov_dict.items():
    ax[0].plot(coordinate_ary, coverage_ary, linewidth=0.25, label = code)
ax[0].set_ylabel('Coverage (x)')
ax[0].set_title('Coverage over chr' + chr + ' for 1x samples')
ax[0].set_xlabel('Genomic position')
ax[0].set_ylim(-0.2, ylim)
ax[0].legend()

for j in range(repeat_mask.shape[0]):
    ax[1].plot(repeat_mask.iloc[j,:2].values, -5*np.ones(2), color='magenta', linewidth=2, alpha=0.3)
ax[1].axis('off')
plt.figtext(0.05, 0.155, 'Repeat Mask')
plt.figtext(.5, 0.1, 
            'Figure 6.' + chr +'. Chr' + chr +' per bin coverage across samples.', ha='center', wrap = True)
plt.savefig('graphs/fig6_per_bin_coverage_chr' + chr +'.png', bbox_inches = "tight", dpi=300)