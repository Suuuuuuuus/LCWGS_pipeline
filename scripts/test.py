import io
import os
import sys
import threading
from queue import Queue
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm
import random
from collections import Counter
import csv
import gzip
import matplotlib.colors as mcolors
from scipy.stats import poisson
import itertools
import collections
import scipy
from scipy.stats import chi2
from scipy.stats import friedmanchisquare
from scipy.stats import studentized_range

def parse_vcf(file, q = None):
    with io.TextIOWrapper(gzip.open(file,'r')) as f:
        lines =[l for l in f if not l.startswith('##')]
        # Identify columns name line and save it into a dict
        # with values as dtype
        dynamic_header_as_key = []
        for liness in f:
            if liness.startswith("#CHROM"):
                dynamic_header_as_key.append(liness)
                # Declare dtypes
        values = [str,int,str,str,str,int,str,str,str,str]
        columns2detype = dict(zip(dynamic_header_as_key,values))

        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype=columns2detype,
            sep='\t'
        ).rename(columns={'#CHROM':'CHROM'})
    if q is None:
        return df
    else:
        q.put(df)

q = Queue()
thread_1 = threading.Thread(target = parse_vcf, args = ('../../GAMCC_oneKG/quilt.chr21.vcf.gz', q))
thread_2 = threading.Thread(target = parse_vcf, args = ('../../GAMCC_oneKG/quilt.chr22.vcf.gz', q))
thread_3 = threading.Thread(target = parse_vcf, args = ('../../GAMCC_oneKG/quilt.chr20.vcf.gz', q))

start = time.time()
thread_1.start()
thread_2.start()
thread_3.start()
thread_1.join()
thread_2.join()
thread_3.join()
end = time.time()
print(end - start)
