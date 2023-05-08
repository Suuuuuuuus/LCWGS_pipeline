import io
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import csv
from scipy.stats import poisson
import itertools
import collections

id = sys.argv[1]
chr = sys.argv[2]
input_path = sys.argv[3]
output_path = sys.argv[4]
bin_size = int(sys.argv[5])

if not os.path.exists(output_path):
    os.makedirs(output_path)

genome_coverage = pd.read_csv(input_path + id + '_bedgraph.txt', header = None, names = ['chr', 'start', 'end', 'cov'], sep = '\t',
                             dtype = {
        'chr': 'string',
        'start': 'Int64',
        'end': 'Int64',
        'cov': 'Int64'
    })

def per_bin_converage(genome_coverage, bin_size = 10000):
    bin_size = bin_size
    start_pos = bin_size/2
    end_pos = genome_coverage['end'].values[-1]-1 - (genome_coverage['end'].values[-1]-1-start_pos)%bin_size

    genome_coverage = genome_coverage[genome_coverage['end'] >= start_pos]
    genome_coverage = genome_coverage[genome_coverage['start'] < end_pos]
    genome_coverage['end'].values[-1] = end_pos-1

    num_bins = int((end_pos - start_pos)/bin_size)
    coordinate_ary = np.arange(bin_size, end_pos, bin_size)
    base_ary = np.zeros(num_bins)
    genome_coverage['fake_end'] = genome_coverage['end'] - start_pos
    genome_coverage['bin_category'] = genome_coverage['fake_end'] // bin_size
    genome_coverage['fake_end'] = genome_coverage['fake_end'].astype(int)
    genome_coverage['bin_category'] = genome_coverage['bin_category'].astype(int)

    prev_end = 0
    prev_bin_category = 0
    for i in range(genome_coverage.shape[0]):
        if genome_coverage['bin_category'].values[i] == prev_bin_category:
            base_ary[prev_bin_category] += genome_coverage['cov'].values[i]*(genome_coverage['fake_end'].values[i] - prev_end)
        elif genome_coverage['bin_category'].values[i] - prev_bin_category == 1:      
            base_ary[prev_bin_category] += genome_coverage['cov'].values[i]*((prev_bin_category+1)*bin_size - prev_end)
            base_ary[prev_bin_category+1] += genome_coverage['cov'].values[i]*(genome_coverage['fake_end'].values[i]- (prev_bin_category+1)*bin_size)
        else:
            base_ary[prev_bin_category] += genome_coverage['cov'].values[i]*((prev_bin_category+1)*bin_size - prev_end)
            base_ary[genome_coverage['bin_category'].values[i]] += genome_coverage['cov'].values[i]*(genome_coverage['fake_end'].values[i]- (genome_coverage['bin_category'].values[i])*bin_size)
            base_ary[prev_bin_category+1:genome_coverage['bin_category'].values[i]] = genome_coverage['cov'].values[i]*bin_size
        prev_bin_category = genome_coverage['bin_category'].values[i]
        prev_end = genome_coverage['fake_end'].values[i]
    base_ary[-1] += genome_coverage['cov'].values[-1]
    return coordinate_ary, base_ary/bin_size

coordinate_ary, base_ary = per_bin_converage(genome_coverage[genome_coverage['chr'] == ('chr'+chr)], bin_size=bin_size)
np.savetxt(output_path + id + '_chr' + chr + '_base.txt', base_ary, newline = ' ', fmt='%1.4f')
np.savetxt(output_path + id + '_chr' + chr + '_coordinate.txt', coordinate_ary, newline = ' ', fmt='%1f')
