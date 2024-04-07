include: "imputation_calculation.smk"
include: "filter_vcf.smk"

include: "auxiliary.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

PANEL_NAME = config['PANEL_NAME']
imp_dir = config['imputation_dir']
case_controls = ['non-malaria_control', 'mild_malaria', 'severe_malaria']
ethnicities = ['fula', 'jola', 'mandinka', 'wollof']
pair = ['lc', 'hc']
axis = ['h', 'v']

rule filter_vcf_all:
    input:
        filtered_vcf = expand(f"results/wip_vcfs/{PANEL_NAME}/high_info/lc.chr{{chr}}.vcf.gz", chr = chromosome)

rule imputation_calculation_all:
    input:
        vcfs_all = expand(imp_dir + 'vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair),
        h_report_all = expand(imp_dir + "impacc/all_samples/by_variant/chr{chr}.h.tsv", chr = chromosome),
        h_impacc_all = expand(imp_dir + "impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome),
        v_report_all = expand(imp_dir + "impacc/all_samples/by_sample/chr{chr}.v.tsv", chr = chromosome),
        v_impacc_all = expand(imp_dir + "impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome),
        r2NRC_h_all = imp_dir + "graphs/all_samples/by_variant/r2_NRC.png",
        ccd_h_all = imp_dir + "graphs/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v_all = imp_dir + "graphs/all_samples/by_sample/r2_NRC.png",
        ccd_v_all = imp_dir + "graphs/all_samples/by_sample/ccd_by_genotype.png",

        vcfs_eth = expand(imp_dir + 'vcf/by_eth/{pair}_vcf/{eth}.{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, eth = ethnicities),
        h_report_eth = expand(imp_dir + "impacc/by_eth/by_variant/{eth}.chr{chr}.h.tsv", eth = ethnicities, chr = chromosome),
        h_impacc_eth = expand(imp_dir + "impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv", eth = ethnicities, chr = chromosome),
        v_report_eth = expand(imp_dir + "impacc/by_eth/by_sample/{eth}.chr{chr}.v.tsv", eth = ethnicities, chr = chromosome),
        v_impacc_eth = expand(imp_dir + "impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv", eth = ethnicities, chr = chromosome),
        r2NRC_h_eth = expand(imp_dir + "graphs/by_eth/by_variant/{eth}.r2_NRC.png", eth = ethnicities),
        ccd_h_eth = expand(imp_dir + "graphs/by_eth/by_variant/{eth}.ccd_by_genotype.png", eth = ethnicities),
        r2NRC_v_eth = expand(imp_dir + "graphs/by_eth/by_sample/{eth}.r2_NRC.png", eth = ethnicities),
        ccd_v_eth = expand(imp_dir + "graphs/by_eth/by_sample/{eth}.ccd_by_genotype.png", eth = ethnicities),

        vcfs_cc = expand(imp_dir + 'vcf/by_cc/{pair}_vcf/{cc}.{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, cc = case_controls),
        h_report_cc = expand(imp_dir + "impacc/by_cc/by_variant/{cc}.chr{chr}.h.tsv", cc = case_controls, chr = chromosome),
        h_impacc_cc = expand(imp_dir + "impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv", cc = case_controls, chr = chromosome),
        v_report_cc = expand(imp_dir + "impacc/by_cc/by_sample/{cc}.chr{chr}.v.tsv", cc = case_controls, chr = chromosome),
        v_impacc_cc = expand(imp_dir + "impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv", cc = case_controls, chr = chromosome),
        r2NRC_h_cc = expand(imp_dir + "graphs/by_cc/by_variant/{cc}.r2_NRC.png", cc = case_controls),
        ccd_h_cc = expand(imp_dir + "graphs/by_cc/by_variant/{cc}.ccd_by_genotype.png", cc = case_controls),
        r2NRC_v_cc = expand(imp_dir + "graphs/by_cc/by_sample/{cc}.r2_NRC.png", cc = case_controls),
        ccd_v_cc = expand(imp_dir + "graphs/by_cc/by_sample/{cc}.ccd_by_genotype.png", cc = case_controls)