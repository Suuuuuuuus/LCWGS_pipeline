configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

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

def get_vcf_concat_lst(region_json, in_prefix):
    REGIONS={}
    for chr in chromosome:
        start = [10000001, 15000001]
        end = [  15000000, 20000000]
        REGIONS[str(chr)]={"start":start, "end":end}

    file = region_json
    if os.path.exists(file):
        with open(file) as json_file:
            REGIONS = json.load(json_file)

    vcfs_to_concat = {}
    vcfs_to_impute = []
    for chr in chromosome:
        start = REGIONS[str(chr)]["start"]
        end = REGIONS[str(chr)]["end"]
        per_chr_vcfs = []
        for i in range(0, start.__len__()):
            regionStart = start[i]
            regionEnd = end[i]
            file = in_prefix + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            per_chr_vcfs.append(file)
            vcfs_to_impute.append(file)
        vcfs_to_concat[str(chr)] = per_chr_vcfs
    return vcfs_to_impute, vcfs_to_concat