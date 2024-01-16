configfile: "pipelines/config.json"

clean_fastq = config['clean_fastq']
reheader = config['reheader']
concatenate = config['concatenate']

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

def read_tsv_as_lst(path): # tsv file should ALWAYS have a single column without header
    return list(pd.read_table(path, header = None, names = ['Code'])['Code'].values)

