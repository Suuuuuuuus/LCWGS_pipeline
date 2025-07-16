import io
import os
import sys
import csv
import gzip
import time
import json
import secrets
import pickle
import multiprocessing
import subprocess
import resource
import itertools
import collections
import sqlite3
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import scipy as sp
import pandas as pd
import statsmodels.api as sm
import random
from collections import Counter
from collections import defaultdict
import copy
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

from scipy.stats import nbinom, norm, geom, beta, poisson
from scipy.special import logsumexp

sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/lcwgsus/')
import lcwgsus
from lcwgsus.variables import *

sys.path.append('/well/band/users/rbx225/GAMCC/scripts/lcSV/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/Low Coverage Data/gamcc/scripts/lcSV/')
from lcSV import *

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

pd.options.mode.chained_assignment = None

def main(regions, sv_df_file, ix, ofile, bin_size = 1000):
    df = pd.read_csv(sv_df_file, sep = '\t')
    chrom, start, L, svtype = df.loc[ix, ['#CHROM', 'POS', 'SVLEN', 'SVTYPE']]
    chromosome = int(chrom.replace('chr', ''))
    start, end = delineate_region(start, L)
    flank = end - start

    cov = load_region_files(regions, chromosome, start, end)
    cov = cov[['position'] + list(cov.columns[cov.columns.str.contains('coverage')])]
    means, variances = normalise_by_flank(cov, start, end, flank)
    samples, coverage = extract_target_cov(cov, start, end)

    results = nonahore(means, variances, coverage, n_recomb = 1000, n_iter = 2000, verbose = False)

    probs, genotypes = results['probs'], results['genotypes']
    info, freq = evaluate_real_model(results)
    haps = results['model_ary'][-1].haps

    outputs = {}
    outputs['chromosome'] = chrom
    outputs['start'] = start
    outputs['end'] = end
    outputs['length'] = L
    outputs['svtype'] = svtype
    outputs['means'] = means
    outputs['variances'] = variances
    outputs['coverage'] = coverage
    outputs['info'] = info
    outputs['freq'] = freq
    outputs['haps'] = haps
    outputs['probs'] = np.round(probs, 4)
    outputs['genotypes'] = genotypes

    with open(ofile, 'wb') as of:
        pickle.dump(outputs, of)
    
if __name__ == "__main__":
    chunk_file = snakemake.params.chunk_file
    with open(chunk_file, 'r') as file:
        regions = json.load(file)

    odir = snakemake.params.odir
    os.makedirs(odir,  exist_ok=True)

    row_ix = int(snakemake.params.row_ix)
    sv_df_file = snakemake.input.sv_df_file
    
    ofile = snakemake.output.pickle

    main(regions, sv_df_file, row_ix, ofile)
