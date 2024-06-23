include: "filter_vcf.smk"
include: "post_gw.smk"
include: "imputation_calculation.smk"

include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

panels = config["panels"][:-1]
PANEL_NAME = config['PANEL_NAME']
imp_dir = config['imputation_dir']
case_controls = ['non-malaria_control', 'mild_malaria', 'severe_malaria']
ethnicities = ['fula', 'jola', 'mandinka', 'wollof']
pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['imputation_dir']

rule post_gw_all:
    input:
        lifted = expand("results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/chr{chr}.dose.vcf.gz", chr = chromosome),
        rejected = expand("results/two-stage-imputation/vanilla/oneKG_hrc/lifted_vcf/chr{chr}.rejected.vcf.gz", chr = chromosome)

rule filter_vcf_all:
    input:
        # lc_chip_site_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        # concat = expand("results/wip_vcfs/{panel_name}/vanilla/chip_sites/lc.vcf.gz", panel_name = panels),
        # PC = "results/wip_vcfs/{panel_name}/vanilla/chip_sites/PCs.eigenvec",
        high_info_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/high_info/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        high_info_high_maf_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        high_info_high_maf_chip_sites_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af_chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        high_info_high_maf_confident_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af_high_conf/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        high_info_high_maf_confident_chip_format = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af_high_conf_chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        # high_info_high_maf_giab_confident_vcf = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af_giab_high_conf/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels),
        # high_info_high_maf_giab_confident_chip_format = expand("results/wip_vcfs/{panel_name}/vanilla/high_info_high_af_giab_high_conf_chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome, panel_name = panels)

rule imputation_calculation_hc_all:
    input:
        vcfs_all = expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', imp_dir = imputation_dir, chr = chromosome, pair = pair),
        h_report_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_report_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_impacc_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        impacc_h = expand("{imp_dir}impacc/all_samples/by_variant/all.h.impacc.tsv", imp_dir = imputation_dir),
        impacc_v = expand("{imp_dir}impacc/all_samples/by_sample/all.v.impacc.tsv", imp_dir = imputation_dir),

        # r2NRC_h_all = expand("{imp_dir}graphs/all_samples/by_variant/r2_NRC.png", imp_dir = imputation_dir),
        # ccd_h_all = expand("{imp_dir}graphs/all_samples/by_variant/ccd_by_genotype.png", imp_dir = imputation_dir),
        # r2NRC_v_all = expand("{imp_dir}graphs/all_samples/by_sample/r2_NRC.png", imp_dir = imputation_dir),
        # ccd_v_all = expand("{imp_dir}graphs/all_samples/by_sample/ccd_by_genotype.png", imp_dir = imputation_dir),

        lc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        hc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        afs = expand("{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv", imp_dir = imputation_dir, chr = chromosome),
        
        tsv = expand("{imp_dir}impacc/all_samples/by_variant/all_r2less0.5.tsv", imp_dir = imputation_dir),
        sumstats = expand("{imp_dir}summary_metrics.tsv", imp_dir = imputation_dir)

rule imputation_calculation_all:
    input:
        vcfs_eth = expand('{imp_dir}vcf/by_eth/{pair}_vcf/{eth}.{pair}.chr{chr}.vcf.gz', imp_dir = imputation_dir, chr = chromosome, pair = pair, eth = ethnicities),
        h_report_eth = expand("{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.tsv", imp_dir = imputation_dir, eth = ethnicities, chr = chromosome),
        h_impacc_eth = expand("{imp_dir}impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, eth = ethnicities, chr = chromosome),
        v_report_eth = expand("{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.tsv", imp_dir = imputation_dir, eth = ethnicities, chr = chromosome),
        v_impacc_eth = expand("{imp_dir}impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, eth = ethnicities, chr = chromosome),
        eth_impacc_h = expand("{imp_dir}impacc/by_cc/by_variant/{eth}.all.h.impacc.tsv", imp_dir = imputation_dir, eth = ethnicities),
        eth_impacc_v = expand("{imp_dir}impacc/by_cc/by_sample/{eth}.all.v.impacc.tsv", imp_dir = imputation_dir, eth = ethnicities),
        # r2NRC_h_eth = expand("{imp_dir}graphs/by_eth/by_variant/{eth}.r2_NRC.png", imp_dir = imputation_dir, eth = ethnicities),
        # ccd_h_eth = expand("{imp_dir}graphs/by_eth/by_variant/{eth}.ccd_by_genotype.png", imp_dir = imputation_dir, eth = ethnicities),
        # r2NRC_v_eth = expand("{imp_dir}graphs/by_eth/by_sample/{eth}.r2_NRC.png", imp_dir = imputation_dir, eth = ethnicities),
        # ccd_v_eth = expand("{imp_dir}graphs/by_eth/by_sample/{eth}.ccd_by_genotype.png", imp_dir = imputation_dir, eth = ethnicities),

        vcfs_cc = expand('{imp_dir}vcf/by_cc/{pair}_vcf/{cc}.{pair}.chr{chr}.vcf.gz', imp_dir = imputation_dir, chr = chromosome, pair = pair, cc = case_controls),
        h_report_cc = expand("{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.tsv", imp_dir = imputation_dir, cc = case_controls, chr = chromosome),
        h_impacc_cc = expand("{imp_dir}impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, cc = case_controls, chr = chromosome),
        v_report_cc = expand("{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.tsv", imp_dir = imputation_dir, cc = case_controls, chr = chromosome),
        v_impacc_cc = expand("{imp_dir}impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, cc = case_controls, chr = chromosome),
        cc_impacc_h = expand("{imp_dir}impacc/by_cc/by_variant/{cc}.all.h.impacc.tsv", imp_dir = imputation_dir, cc = case_controls),
        cc_impacc_v = expand("{imp_dir}impacc/by_cc/by_sample/{cc}.all.v.impacc.tsv", imp_dir = imputation_dir, cc = case_controls),
        # r2NRC_h_cc = expand("{imp_dir}graphs/by_cc/by_variant/{cc}.r2_NRC.png", imp_dir = imputation_dir, cc = case_controls),
        # ccd_h_cc = expand("{imp_dir}graphs/by_cc/by_variant/{cc}.ccd_by_genotype.png", imp_dir = imputation_dir, cc = case_controls),
        # r2NRC_v_cc = expand("{imp_dir}graphs/by_cc/by_sample/{cc}.r2_NRC.png", imp_dir = imputation_dir, cc = case_controls),
        # ccd_v_cc = expand("{imp_dir}graphs/by_cc/by_sample/{cc}.ccd_by_genotype.png", imp_dir = imputation_dir, cc = case_controls),
        # vcfs_all = expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', imp_dir = imputation_dir, chr = chromosome, pair = pair),
        h_report_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv", imp_dir = imputation_dir, chr = chromosome),
        h_impacc_all = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_report_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv", imp_dir = imputation_dir, chr = chromosome),
        v_impacc_all = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", imp_dir = imputation_dir, chr = chromosome),
        impacc_h = expand("{imp_dir}impacc/by_cc/by_variant/all.h.impacc.tsv", imp_dir = imputation_dir),
        impacc_v = expand("{imp_dir}impacc/by_cc/by_sample/all.v.impacc.tsv", imp_dir = imputation_dir),
        # r2NRC_h_all = expand("{imp_dir}graphs/all_samples/by_variant/r2_NRC.png", imp_dir = imputation_dir),
        # ccd_h_all = expand("{imp_dir}graphs/all_samples/by_variant/ccd_by_genotype.png", imp_dir = imputation_dir),
        # r2NRC_v_all = expand("{imp_dir}graphs/all_samples/by_sample/r2_NRC.png", imp_dir = imputation_dir),
        # ccd_v_all = expand("{imp_dir}graphs/all_samples/by_sample/ccd_by_genotype.png", imp_dir = imputation_dir),

        # lc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        # hc_vcfs = expand("{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz", imp_dir = imputation_dir, chr = chromosome),
        # afs = expand("{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv", imp_dir = imputation_dir, chr = chromosome),
        
        sumstats = expand("{imp_dir}summary_metrics_all.tsv", imp_dir = imputation_dir)