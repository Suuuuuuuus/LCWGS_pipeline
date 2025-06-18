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

def main(odir):
    N = 210
    nb_vars = [100, 500, 1000, 2000, 5000]
    Ls = [5, 10, 20, 30] # translate to 1k, 3k, 6.6k, 10k SVLEN
    fs = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    binsize = 1000

    mean_coverage = 1.21
    sd_coverage = 0.12

    n_recomb = 200
    n_iter = 500
    verbose = False

    eval_dict = {}

    for v in nb_vars:
        for l in Ls:
            for f1 in fs:  
                h1 = np.ones(l)
                h2 = h1.copy()
                h2[int(l/3):int(2*l/3)] = 0
                model = SVModel([h1,h2], [1-f1,f1])
                training, coverage, true_gt = simulate_coverotron(model, N, l, mean_coverage, sd_coverage, v)

                np.savetxt(f'{odir}{v}-{l}-{f1}-training.txt', training, delimiter='\t', fmt = "%d")
                np.savetxt(f'{odir}{v}-{l}-{f1}-coverage.txt', coverage, delimiter='\t', fmt = "%d")
                np.savetxt(f'{odir}{v}-{l}-{f1}-truth.txt', true_gt, delimiter='\t', fmt = "%d")

                means = np.mean(training, axis = 0)
                variances = np.var(training, axis = 0, ddof = 1)
                result_dict = nonahore(means, variances, coverage, n_recomb = n_recomb, n_iter = n_iter, verbose = verbose)

                concordance, info, freq = evaluate_sim_model(result_dict, h2, true_gt)
                eval_dict[f'{v}-{l}-{f1}'] = [concordance, info, freq]

    pickle_ofile = f'{odir}eval.pickle'
    with open(pickle_ofile, 'wb') as of:
        pickle.dump(eval_dict, of)
    
if __name__ == "__main__":
    odir = snakemake.params.odir
    os.makedirs(odir,  exist_ok=True)
    
    main(odir)
