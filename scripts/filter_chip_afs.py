import io
import os
import sys
import csv
import gzip
import time
import secrets
import multiprocessing
import resource
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import random
from collections import Counter
import seaborn as sns
import matplotlib.colors as mcolors
import itertools
import collections
sys.path.append('.')
import lcwgSus

from scipy.stats import poisson
from scipy.stats import chi2
from scipy.stats import friedmanchisquare
from scipy.stats import studentized_range
pd.options.mode.chained_assignment = None

chip_vcf = "../results/chip/filtered_snps.vcf.gz"
df = lcwgSus.read_vcf(chip_vcf)
colnames = lcwgSus.read_metadata(chip_vcf)

chromosomes = [i for i in range(1,23)]
mafs = ["/well/band/users/rbx225/recyclable_files/AFs/ggvp_AFs/ggvp_AF_chr" + str(i) + ".txt" for i in chromosomes]
af = lcwgSus.multi_read_af(chromosomes, mafs)
af = af.rename(columns = {'MAF': 'prop'})
chip_af = lcwgSus.calculate_af(df)

def filter_afs(df1, df2, diff = 0.2, z_score = None):
    # df1 is the main vcf in which afs are to be filtered out
    # df2 is the ref panel afs
    # Either filter by z-score (suggested 2 sds so 1.96 or diff=0.2)
    res = pd.merge(df1, df2, on = ['chr', 'pos', 'ref', 'alt'])
    if z_score is not None:
        res = res[(res['prop_y'] != 0) & (res['prop_y'] != 1)]
        res['z'] = (res['prop_x'] - res['prop_y'])/np.sqrt(res['prop_y']*(1-res['prop_y']))
        res = res[abs(res['z']) <= z_score]
        return res.drop(columns = ['prop_x', 'prop_y', 'z'])
    else:
        res = res[abs(res['prop_x'] - res['prop_y']) < diff]
        return res.drop(columns = ['prop_y']).rename(columns = {'prop_x': 'prop'})

chip_res = filter_afs(chip_af, af)
lcwgSus.save_vcf(chip_res, colnames, '../results/chip/filtered_snps.qced.vcf.gz')
lcwgSus.plot_afs(chip_res, af, save_fig = True, outdir = '../graphs/')
