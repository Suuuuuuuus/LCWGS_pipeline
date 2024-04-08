import importlib
import io
import os
import sys
import csv
import gzip
import time
import secrets
import multiprocessing
import subprocess
import resource
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import random
from collections import Counter
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter
import itertools
import collections
sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
import lcwgsus

from scipy.stats import poisson
from scipy.stats import chi2
from scipy.stats import friedmanchisquare
from scipy.stats import studentized_range
pd.options.mode.chained_assignment = None

chromosome = [str(i) for i in range(1,23)]
common_cols = ['chr', 'pos', 'ref', 'alt']

imp_dir = "/well/band/users/rbx225/GAMCC/results/imputation_comparison/oneKG/lc_chip_typed_high_info/"

for c in chromosome:
    quilt_vcf = imp_dir + "vcf/all_samples/lc_vcf/lc.chr" + c + ".vcf.gz"
    chip_vcf = imp_dir + "vcf/all_samples/hc_vcf/hc.chr" + c + ".vcf.gz"
    af_txt = "/well/band/users/rbx225/GAMCC/data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr" + c + ".txt"

    common_savedir = imp_dir + "vcf/all_samples/"
    common_outdir = imp_dir + "impacc/all_samples/"

    chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt, 
                                                             save_vcfs = True, 
                                                             lc_vcf_outdir = common_savedir + "filtered_vcfs/", 
                                                             hc_vcf_outdir = common_savedir + "filtered_vcfs/", 
                                                             af_outdir = common_savedir + "af/",
                                                             lc_vcf_name = "lc.chr" +c+ ".vcf.gz", 
                                                             hc_vcf_name = "hc.chr" +c+ ".vcf.gz", 
                                                             af_name = "af.chr" +c+ ".tsv")