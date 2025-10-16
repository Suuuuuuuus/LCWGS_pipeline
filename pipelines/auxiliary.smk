configfile: "pipelines/config.json"
home_dir = config['home_dir']

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append(f"{home_dir}software/lcwgsus/")
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

# This utility extracts all values in a dictionary and combine them in a flattened list
def convert_dict_to_lst(dictionary):
    return [item for sublist in dictionary.values() for item in sublist]

def get_vcf_concat_lst(region_json, ref_prefix, vcf_prefix, suffix = '.vcf.gz'):
    REGIONS = {}
    for chr in chromosome:
        start = [10000001, 15000001]
        end = [  15000000, 20000000]
        REGIONS[str(chr)] = {"start": start, "end": end}

    if os.path.exists(region_json):
        with open(region_json) as json_file:
            REGIONS = json.load(json_file)

    regions_to_prep = []
    vcfs_to_concat = {}
    for chr in chromosome:
        start = REGIONS[str(chr)]["start"]
        end = REGIONS[str(chr)]["end"]
        per_chr_vcfs = []
        for i in range(0, start.__len__()):
            regionStart = start[i]
            regionEnd = end[i]
            file = ref_prefix + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
            regions_to_prep.append(file)
            file = vcf_prefix + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + suffix
            per_chr_vcfs.append(file)
        vcfs_to_concat[str(chr)] = per_chr_vcfs

    vcfs_to_impute = convert_dict_to_lst(vcfs_to_concat)
    return regions_to_prep, vcfs_to_impute, vcfs_to_concat

def get_bam_concat_lst(region_json, sample_list, ref_prefix, vcf_prefix, bam_suffix = '.bam', output_suffix = '.tsv.gz'):
    REGIONS = {}
    for chr in chromosome:
        start = [10000001, 15000001]
        end = [  15000000, 20000000]
        REGIONS[str(chr)] = {"start": start, "end": end}

    if os.path.exists(region_json):
        with open(region_json) as json_file:
            REGIONS = json.load(json_file)

    regions_to_prep = {}
    cov_files_all = []
    for chr in chromosome:
        start = REGIONS[str(chr)]["start"]
        end = REGIONS[str(chr)]["end"]
        for i in range(0, start.__len__()):
            regionStart = start[i]
            regionEnd = end[i]
            bamlist = []
            for s in sample_list:
                file = ref_prefix + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + '/' + s + bam_suffix
                bamlist.append(file)
            regions_to_prep[f'{chr}.{regionStart}.{regionEnd}'] = bamlist
            file = vcf_prefix + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + output_suffix
            cov_files_all.append(file)
            
    return regions_to_prep, cov_files_all
