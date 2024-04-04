include: "auxiliary.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_chip = read_tsv_as_lst(config['samples_chip'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
chromosome = [i for i in range(1,23)]

imp_dir = config['imputation_dir']
case_controls = ['non-maleria_control', 'mild_malaria', 'severe_malaria']
ethnicities = ['fula', 'jola', 'mandinka', 'wollof']
pair = ['lc', 'hc']
axis = ['h', 'v']

rule preprocess_vcf_in_working_dir:
    input:
        lc_vcf_dir = config['lc_vcf_dir'],
        hc_vcf_dir = config['hc_vcf_dir']
    output:
        vcfs = expand(imp_dir + 'vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair),
        cc_vcfs = expand(imp_dir + 'vcf/by_cc/{pair}_vcf/{cc}.{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, cc = case_controls),
        eth_vcfs = expand(imp_dir + 'vcf/by_eth/{pair}_vcf/{eth}.{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, eth = ethnicities)
    # params:
        # mild_malaria_lc = config['mild_malaria_lc']
        # non_malaria_control_lc = config['non_malaria_control_lc']
        # severe_malaria_lc = config['severe_malaria_lc']
        # mild_malaria_chip = config['mild_malaria_chip']
        # non_malaria_control_chip = config['non_malaria_control_chip']
        # severe_malaria_chip = config['severe_malaria_chip']

        # fula_lc = config['fula_lc']
        # jola_lc = config['jola_lc']
        # mandinka_lc = config['mandinka_lc']
        # wollof_lc = config['wollof_lc']
        # fula_chip = config['fula_chip']
        # jola_chip = config['jola_chip']
        # mandinka_chip = config['mandinka_chip']
        # wollof_chip = config['wollof_chip']
    shell: """
        mkdir -p "{imp_dir}vcf/" "{imp_dir}impacc/" "{imp_dir}graphs/"
        mkdir -p "{imp_dir}vcf/all_samples/lc_vcf/" "{imp_dir}vcf/all_samples/hc_vcf/"
        mkdir -p "{imp_dir}vcf/by_cc/lc_vcf/" "{imp_dir}vcf/by_cc/hc_vcf/"
        mkdir -p "{imp_dir}vcf/by_eth/lc_vcf/" "{imp_dir}vcf/by_eth/hc_vcf/"

        declare -a eth=("jola" "fula" "mandinka" "wollof")
        declare -a cc=("mild_malaria" "non-malaria_control" "severe_malaria")
        declare -a pair=("lc" "hc")
        
        for i in {1..22}
        do 
            cp $"{input.lc_vcf_dir}*chr$i*.gz" "{imp_dir}vcf/all_samples/lc_vcf/lc.chr$i.vcf.gz"
            cp $"{input.hc_vcf_dir}*chr$i*.gz" "{imp_dir}vcf/all_samples/hc_vcf/hc.chr$i.vcf.gz"
            
            for e in "${{eth[@]}}"
            do
                for p in "${{pair[@]}}"
                do
                    bcftools view -S "file_lsts/samples_subset/by_ethnicity/"$e"_samples_"$p".tsv" -Oz -o "{imp_dir}vcf/by_eth/"$p"_vcf/"$e"."$p".chr"$i".vcf.gz" "{imp_dir}vcf/all_samples/"$p"_vcf/"$p".chr"$i".vcf.gz"
                done
            done
            
            for c in "${{cc[@]}}"
            do
                for p in "${{pair[@]}}"
                do
                    bcftools view -S "file_lsts/samples_subset/by_case_control/"$c"_samples_"$p".tsv" -Oz -o "{imp_dir}vcf/by_cc/"$p"_vcf/"$c"."$p".chr"$i".vcf.gz" "{imp_dir}vcf/all_samples/"$p"_vcf/"$p".chr"$i".vcf.gz"
                done
            done

        done
    """

rule calculate_imputation_accuracy_all:
    input:
        quilt_vcf = imp_dir + 'vcf/all_samples/lc_vcf/lc.chr{chr}.vcf.gz',
        chip_vcf = imp_dir + 'vcf/all_samples/hc_vcf/hc.chr{chr}.vcf.gz',
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = imp_dir + "impacc/all_samples/by_variant/chr{chr}.h.tsv",
        h_impacc = imp_dir + "impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv",
        v_report = imp_dir + "impacc/all_samples/by_sample/chr{chr}.v.tsv",
        v_impacc = imp_dir + "impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker'],
        common_outdir = imp_dir + "impacc/all_samples/"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = 'chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = 'chr' + wildcards.chr +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy_all:
    input:
        h_impaccs = expand(imp_dir + "impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome),
        v_impaccs = expand(imp_dir + "impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome)
    output:
        r2NRC_h = imp_dir + "graphs/all_samples/by_variant/r2_NRC.png",
        ccd_h = imp_dir + "graphs/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = imp_dir + "graphs/all_samples/by_sample/r2_NRC.png",
        ccd_v = imp_dir + "graphs/all_samples/by_sample/ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = imp_dir + "graphs/all_samples/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.common_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = "r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = "ccd_by_genotype.png")

        v_lst = input.v_impaccs
        outdir_v = params.common_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = "r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = "ccd_by_genotype.png")

rule calculate_imputation_accuracy_by_eth:
    input:
        quilt_vcf = imp_dir + 'vcf/by_eth/lc_vcf/{eth}.lc.chr{chr}.vcf.gz',
        chip_vcf = imp_dir + 'vcf/by_eth/hc_vcf/{eth}.hc.chr{chr}.vcf.gz',
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = imp_dir + "impacc/by_eth/by_variant/{eth}.chr{chr}.h.tsv",
        h_impacc = imp_dir + "impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv",
        v_report = imp_dir + "impacc/by_eth/by_sample/{eth}.chr{chr}.v.tsv",
        v_impacc = imp_dir + "impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker'],
        common_outdir = imp_dir + "impacc/by_eth/"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = wildcards.eth + '.chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.eth + '.chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy_by_eth:
    input:
        h_impaccs = expand(imp_dir + "impacc/by_eth/by_variant/{eth}.chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand(imp_dir + "impacc/by_eth/by_sample/{eth}.chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        r2NRC_h = imp_dir + "graphs/by_eth/by_variant/{eth}.r2_NRC.png",
        ccd_h = imp_dir + "graphs/by_eth/by_variant/{eth}.ccd_by_genotype.png",
        r2NRC_v = imp_dir + "graphs/by_eth/by_sample/{eth}.r2_NRC.png",
        ccd_v = imp_dir + "graphs/by_eth/by_sample/{eth}.ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = imp_dir + "graphs/by_eth/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.common_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = wildcards.eth + ".r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = wildcards.eth + ".ccd_by_genotype.png")

        v_lst = input.v_impaccs
        outdir_v = params.common_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = wildcards.eth + ".r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = wildcards.eth + ".ccd_by_genotype.png")

rule calculate_imputation_accuracy_by_cc:
    input:
        quilt_vcf = imp_dir + 'vcf/by_cc/lc_vcf/{cc}.lc.chr{chr}.vcf.gz',
        chip_vcf = imp_dir + 'vcf/by_cc/hc_vcf/{cc}.hc.chr{chr}.vcf.gz',
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = imp_dir + "impacc/by_cc/by_variant/{cc}.chr{chr}.h.tsv",
        h_impacc = imp_dir + "impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv",
        v_report = imp_dir + "impacc/by_cc/by_sample/{cc}.chr{chr}.v.tsv",
        v_impacc = imp_dir + "impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker'],
        common_outdir = imp_dir + "impacc/by_cc/"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = wildcards.cc + '.chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = wildcards.cc + '.chr' + wildcards.chr +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy_by_cc:
    input:
        h_impaccs = expand(imp_dir + "impacc/by_cc/by_variant/{cc}.chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand(imp_dir + "impacc/by_cc/by_sample/{cc}.chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        r2NRC_h = imp_dir + "graphs/by_cc/by_variant/{cc}.r2_NRC.png",
        ccd_h = imp_dir + "graphs/by_cc/by_variant/{cc}.ccd_by_genotype.png",
        r2NRC_v = imp_dir + "graphs/by_cc/by_sample/{cc}.r2_NRC.png",
        ccd_v = imp_dir + "graphs/by_cc/by_sample/{cc}.ccd_by_genotype.png"
    resources:
        mem = '30G'
    params:
        common_outdir = imp_dir + "graphs/by_cc/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.common_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = wildcards.cc + ".r2_NRC.png")
        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = wildcards.cc + ".ccd_by_genotype.png")

        v_lst = input.v_impaccs
        outdir_v = params.common_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = wildcards.cc + ".r2_NRC.png")
        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = wildcards.cc + ".ccd_by_genotype.png")
