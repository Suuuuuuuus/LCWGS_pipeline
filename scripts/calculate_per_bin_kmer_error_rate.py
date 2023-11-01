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

code = sys.argv[1]
strand = sys.argv[2]
file = sys.argv[3]

def calculate_per_bin_kmer_error_rate(code, strand, file, intervals = [50,100,150,200,250,300,350,400,450,500,550,600]):
    result = pd.read_csv(file, sep = '\t', usecols = ['read_id', 'number_of_kmers_at_threshold', 'number_of_solid_kmers_at_threshold'])
    result['read_id'] = result['read_id'].str.split().str[0]
    bins = pd.read_csv('results/kmer/'+code+'/fragment_length.tsv', sep = '\t', header = None, names = ['read_id', 'bin'], dtype={'read_id': 'string', 'bin': 'int'})
    df = pd.merge(result, bins, on='read_id', how="left")

    intervals = intervals
    error_rate_ary=[]
    for i in intervals:
        error_rate = df.loc[df['bin'] == i, 'number_of_solid_kmers_at_threshold'].sum()/df.loc[df['bin'] == i, 'number_of_kmers_at_threshold'].sum()
        error_rate_ary.append(error_rate)
    error_rate = df.loc[df['bin']>intervals[-1], 'number_of_solid_kmers_at_threshold'].sum()/df.loc[df['bin']>intervals[-1], 'number_of_kmers_at_threshold'].sum()
    error_rate_ary.append(error_rate)

    res = pd.DataFrame({code: error_rate_ary}).T
    res.to_csv('results/kmer/'+code+'/read'+strand+'/per_bin_kmer_error_rate_read'+strand+'.txt', header=False, sep='\t')

calculate_per_bin_kmer_error_rate(code, strand, file)
