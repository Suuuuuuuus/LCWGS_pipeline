configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

def read_tsv_as_lst(path): # tsv file should ALWAYS have a single column without header
    if os.path.exists(path):
        return list(pd.read_table(path, header = None, names = ['Code'])['Code'].values)
    else:
        return []

def read_tsv_as_dict(samples, split_path_prefix, split_path_postfix):
    # By default `samples` take a list of samples; if `samples` comes as a path of a tsv file, this could also be detected and dealt accordingly
    if type(samples) == str:
        samples = read_tsv_as_lst(samples)
    dict = {}
    for i in samples:
        path = split_path_prefix + i + split_path_postfix
        dict[i] = read_tsv_as_lst(path)
    return dict