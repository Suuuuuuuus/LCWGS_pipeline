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
chip_df = lcwgSus.read_vcf(chip_vcf)
colnames = lcwgSus.read_metadata(chip_vcf)

chromosomes = [i for i in range(1,23)]
mafs = ["/well/band/users/rbx225/recyclable_files/AFs/ggvp_AFs/ggvp_AF_chr" + str(i) + ".txt" for i in chromosomes]
af = lcwgSus.multi_read_af(chromosomes, mafs)
af = af.rename(columns = {'MAF': 'prop'})
chip_df = lcwgSus.calculate_af(chip_df, drop = False)

chip_res = lcwgSus.filter_afs(chip_df, af)
#lcwgSus.plot_afs(chip_res[['chr', 'pos', 'ref', 'alt', 'prop']], af, save_fig = True, outdir = '../graphs/')
chip_res = chip_res.drop(columns = 'prop')
lcwgSus.save_vcf(chip_res, colnames, '../results/chip/filtered_snps.qced.vcf.gz')



