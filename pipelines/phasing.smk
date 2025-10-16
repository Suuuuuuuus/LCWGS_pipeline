configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
sys.path.append(f'{home_dir}software/QUILT_sus/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *
from hla_phase import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

'''

rule phase_GAMCC_alleles:
    input:
        phased_vcf_file = "results/imputation/vcfs/malariaGen_v1_b38/quilt.chr6.vcf.gz"
    output:
        # html = "results/phasing/html/GAMCC-{gene}.html",
        phase_df = "results/phasing/phased_dfs/GAMCC-{gene}.tsv"
    resources: mem = '30G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv',
        reference_allele_file = '/well/band/users/rbx225/recyclable_files/hla/b38_reference_alleles.tsv',
        gm_tsv_names = 'data/sample_tsvs/fv_gm_names.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')

        gamcc_hla = lcwgsus.read_hla_direct_sequencing(retain = 'fv', unique_two_field = False)
        gamcc_hla = gamcc_hla[['SampleID', 'Locus', 'Two field1', 'Two field2']].reset_index(drop = True)

        colnames = ['Sample ID'] + [label for g in HLA_GENES for label in [f'HLA-{g} 1', f'HLA-{g} 2']]
        hlatypes = pd.DataFrame(columns = colnames)
        for s in gamcc_hla['SampleID'].unique():
            tmp = gamcc_hla[gamcc_hla['SampleID'] == s]
            row = [s] + tmp[['Two field1', 'Two field2']].values.ravel().tolist()
            hlatypes.loc[len(hlatypes)] = row
    
        reference_allele_ary = np.array(lcwgsus.read_tsv_as_lst(params.reference_allele_file))

        subset_vcf_samples = lcwgsus.read_tsv_as_lst(params.gm_tsv_names)
        subset_vcf_samples = ','.join(subset_vcf_samples)

        sample_linker = pd.read_csv(SAMPLE_LINKER_FILE)
        sample_linker = {k:v for k, v in zip(sample_linker['Sample_Name'], sample_linker['Chip_Name'])}

        return_dict = phase_hla_on_haplotypes(gene = wildcards.gene, 
                            ipd_gen_file_dir = params.ipd_gen_file_dir, 
                            hla_gene_information = hla_gene_information, 
                            hlatypes = hlatypes,
                            phased_vcf = input.phased_vcf_file, 
                            reference_allele_ary = reference_allele_ary, 
                            strict_snp_filter = False,
                            read_from_QUILT = True, 
                            subset_vcf_samples = subset_vcf_samples,
                            sample_linker = sample_linker)
        # hlatypes = return_dict['hlatypes']

        # individual = 'GAM951638'
        # ix = hlatypes.index[hlatypes['Sample ID'] == individual][0]
        # display_indices = np.arange(10)

        # res = visualise_phase(wildcards.gene, ix, hlatypes, return_dict, both_het = True)
        # compare_phase(display_indices, res, save_html = True, save_name = output.html)
        df = return_dict['phase_df'][['Sample', 'allele1', 'allele2']]
        df.columns = ['Sample ID', f'HLA-{wildcards.gene} 1', f'HLA-{wildcards.gene} 2']
        df.to_csv(output.phase_df, sep = '\t', index = False, header = True)

vcf_versions = ['30x', 'phase3_b38']

rule phase_1KG_alleles:
    input:
        samples_file = "results/hla/imputation/ref_panel/auxiliary_files/oneKG.samples",
        hlatypes_file = "results/hla/imputation/ref_panel/auxiliary_files/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        phased_vcf_file = "/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_{vcf_versions}/oneKG.chr6.vcf.gz"
    output:
        # html = "results/phasing/html/oneKG-{filter}-{gene}.html",
        phase_df = "results/phasing/phased_dfs/oneKG_{vcf_versions}-{filter}-{gene}.tsv"
    resources:
        mem = '30G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv',
        reference_allele_file = '/well/band/users/rbx225/recyclable_files/hla/b38_reference_alleles.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        hlatypes = pd.read_csv(input.hlatypes_file, sep = '\t')

        ref_samples = pd.read_csv(input.samples_file, sep = ' ')
        ref_samples_removed = ref_samples[~ref_samples['SAMPLE'].isin(hlatypes['Sample ID'].tolist())]
        samples_to_remove = ref_samples_removed['SAMPLE'].tolist()
        hlatypes = hlatypes[~hlatypes['Sample ID'].isin(samples_to_remove)].sort_values(by = 'Sample ID').reset_index(drop = True)

        reference_allele_ary = np.array(lcwgsus.read_tsv_as_lst(params.reference_allele_file))
        strict_snp_filter = True if wildcards.filter == 'strict' else False

        return_dict = phase_hla_on_haplotypes(gene = wildcards.gene, 
                                    ipd_gen_file_dir = params.ipd_gen_file_dir, 
                                    hla_gene_information = hla_gene_information,
                                    hlatypes = hlatypes,
                                    phased_vcf = input.phased_vcf_file, 
                                    reference_allele_ary = reference_allele_ary, 
                                    strict_snp_filter = strict_snp_filter,
                                    read_from_QUILT = False, 
                                    subset_vcf_samples = None,
                                    sample_linker = None)
        # hlatypes = return_dict['hlatypes']

        # individual = 'NA12878'
        # ix = hlatypes.index[hlatypes['Sample ID'] == individual][0]
        # display_indices = np.arange(10)

        # res = visualise_phase(wildcards.gene, ix, hlatypes, return_dict, both_het = True)
        # compare_phase(display_indices, res, save_html = True, save_name = output.html)
        df = return_dict['phase_df'][['Sample', 'allele1', 'allele2']]
        df.columns = ['Sample ID', f'HLA-{wildcards.gene} 1', f'HLA-{wildcards.gene} 2']
        df.to_csv(output.phase_df, sep = '\t', index = False, header = True)

rule prepare_1kg_HLA_vcf:
    input:
        hlatypes = "/well/band/users/rbx225/GAMCC/results/hla/imputation/ref_panel/auxiliary_files/HLA_types_2568_for_phasing.txt",
        phased_vcf = "/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/oneKG.chr6.vcf.gz"
    output:
        vcf = "results/phasing/HLA_1KG_BEAGLE/unphased.1KG.chr6.vcf.gz"
    resources: mem = '80G'
    threads: 8
    shell: """
        python scripts/prepare_1KG_BEAGLE_phasing.py
    """

rule prepare_GAMCC_HLA_vcf:
    input:
        phased_vcf = "results/imputation/vcfs/malariaGen_v1_b38/quilt.chr6.vcf.gz"
    output:
        vcf = "results/phasing/HLA_GAMCC_BEAGLE/unphased.GAMCC.chr6.vcf.gz"
    resources: mem = '40G'
    threads: 4
    shell: """
        python scripts/prepare_GAMCC_BEAGLE_phasing.py
    """
'''
rule beagle_phasing:
    input:
        vcf = "results/phasing/HLA_{study}_BEAGLE/unphased.{study}.chr6.vcf.gz"
    output:
        phased_vcf = "results/phasing/HLA_{study}_BEAGLE/phased.{study}.chr6.vcf.gz",
        tmp_vcf = temp("results/phasing/HLA_{study}_BEAGLE/tmp.phased.{study}.chr6.vcf.gz")
    params:
        beagle = tools['beagle'],
        recomb_map = "/well/band/users/rbx225/recyclable_files/plink_recomb_maps_b38/chr6.map",
        output_prefix = "results/phasing/HLA_{study}_BEAGLE/phased.{study}.chr6"
    resources: mem = '40G'
    threads: 4
    shell: """
        mkdir -p results/phasing/HLA_{wildcards.study}_BEAGLE/
        {params.beagle} gt={input.vcf} map={params.recomb_map} out={params.output_prefix}

        mv {output.phased_vcf} {output.tmp_vcf}
        tabix -f {output.tmp_vcf}
        bcftools view -Oz -o {output.phased_vcf} {output.tmp_vcf}
        tabix -f {output.phased_vcf}
    """

rule extract_beagle_phase_vcf:
    input:
        vcf = "results/phasing/HLA_{study}_BEAGLE/phased.{study}.chr6.vcf.gz"
    output:
        vcf = temp("results/phasing/HLA_{study}_BEAGLE/tmp/phased.{gene}.{study}.chr6.vcf.gz")
    resources: mem = '40G'
    threads: 4
    shell: """ 
        mkdir -p results/phasing/HLA_{wildcards.study}_BEAGLE/tmp/
        
        bcftools filter -i 'ID ~ "HLA_{wildcards.gene}"' -Oz -o {output.vcf} {input.vcf}
    """

rule extract_beagle_phase_vcf_hlatypes:
    input:
        vcf = "results/phasing/HLA_{study}_BEAGLE/tmp/phased.{gene}.{study}.chr6.vcf.gz"
    output:
        tsv = temp("results/phasing/HLA_{study}_BEAGLE/tmp/{gene}.semiphased.{study}.tsv")
    resources: mem = '80G'
    threads: 8
    run: 
        vcf = lcwgsus.read_vcf(input.vcf)
        samples = vcf.columns[9:].tolist()
        vcf = vcf.drop(columns = ['chr', 'pos', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT']).reset_index(drop = True).T

        beagle_hla = pd.DataFrame({'Sample ID': samples})
        beagle_hla[f'HLA-{wildcards.gene} 1'] = 'N/A'
        beagle_hla[f'HLA-{wildcards.gene} 2'] = 'N/A'

        for i,s in enumerate(samples):
            vcf_row = vcf.loc[s, :]
            alleles_idx = vcf_row.index[vcf_row.str.contains('1')]
            if len(alleles_idx) != 0:
                if (len(alleles_idx) == 2) and (vcf.loc[s,alleles_idx[0]] == vcf.loc[s,alleles_idx[1]]): # Means beagle phase the two HLA alleles to one haplotype
                    pass
                else:
                    for j in alleles_idx:
                        name = vcf.loc['ID',j]
                        genotype = vcf.loc[s,j]
                        twofield = name.split('*')[1]
                        if genotype == '1|0':
                            beagle_hla.loc[i, f'HLA-{wildcards.gene} 1'] = twofield
                        elif genotype == '0|1':
                            beagle_hla.loc[i, f'HLA-{wildcards.gene} 2'] = twofield
                        elif genotype == '1|1':
                            beagle_hla.loc[i, f'HLA-{wildcards.gene} 1'] = twofield
                            beagle_hla.loc[i, f'HLA-{wildcards.gene} 2'] = twofield
                        else:
                            pass
        beagle_hla.to_csv(output.tsv, header = True, index = False, sep = '\t')

rule adjust_beagle_phasing:
    input:
        unphased_vcf = "results/phasing/HLA_{study}_BEAGLE/unphased.{study}.chr6.vcf.gz",
        phased_vcf = "results/phasing/HLA_{study}_BEAGLE/phased.{study}.chr6.vcf.gz",
        beagle_semiphased = "results/phasing/HLA_{study}_BEAGLE/tmp/{gene}.semiphased.{study}.tsv"
    output:
        beagle_phased = "results/phasing/HLA_{study}_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.{study}.tsv"
    params:
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv'
    resources: mem = '40G'
    threads: 4
    run: 
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        beagle_hla = pd.read_csv(input.beagle_semiphased, sep = '\t')
        
        gene = wildcards.gene
        s = wildcards.sample
        beagle_hla = beagle_hla[beagle_hla['Sample ID'] == s].reset_index(drop = True)

        if pd.isna(beagle_hla.iloc[0, 1]) and pd.isna(beagle_hla.iloc[0, 2]):
            beagle_hla.iloc[0, 1] = 'N/A'
            beagle_hla.iloc[0, 2] = 'N/A'
            beagle_hla.to_csv(output.beagle_phased, header = True, index = False, sep = '\t')
        elif (beagle_hla.iloc[0, 1] == 'N/A') and (beagle_hla.iloc[0, 2] == 'N/A'):
            beagle_hla.iloc[0, 1] = 'N/A'
            beagle_hla.iloc[0, 2] = 'N/A'
            beagle_hla.to_csv(output.beagle_phased, header = True, index = False, sep = '\t')
        elif beagle_hla.iloc[0, 1] == beagle_hla.iloc[0, 2]:
            beagle_hla.iloc[0, 1] = 'N/A'
            beagle_hla.iloc[0, 2] = 'N/A'
            beagle_hla.to_csv(output.beagle_phased, header = True, index = False, sep = '\t')
        else:
            start = hla_gene_information[hla_gene_information['Name'] == f'HLA-{gene}']['Start'].values[0]
            end = hla_gene_information[hla_gene_information['Name'] == f'HLA-{gene}']['End'].values[0]

            idx = []
            multiplier = 1
            while (len(idx) == 0) and (multiplier < 10):
                tmp_new_up = read_vcf(start = max(25000000, start - multiplier*10000), end = max(25000000, start - (multiplier - 1)*10000), phased_vcf = input.phased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                tmp_old_up = read_vcf(start = max(25000000, start - multiplier*10000), end = max(25000000, start - (multiplier - 1)*10000), phased_vcf = input.unphased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                tmp_new_up = pd.merge(tmp_new_up, tmp_old_up[['snp', 'pos']], how = 'inner')
                tmp_old_up = tmp_old_up.reset_index(drop = True)
                
                idx_lst = tmp_new_up.index[tmp_new_up[s].isin(['0|1', '1|0'])]
                if len(idx_lst) != 0:
                    idx.append(idx_lst[-1])
                    new = tmp_new_up
                    old = tmp_old_up
                else:
                    tmp_new_down = read_vcf(start = min(34000000, end + (multiplier - 1)*10000), end = min(34000000, end + multiplier*10000), phased_vcf = input.phased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                    tmp_old_down = read_vcf(start = min(34000000, end + (multiplier - 1)*10000), end = min(34000000, end + multiplier*10000), phased_vcf = input.unphased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                    tmp_new_down = pd.merge(tmp_new_down, tmp_old_down[['snp', 'pos']], how = 'inner')
                    tmp_old_down = tmp_old_down.reset_index(drop = True)

                    idx_lst = tmp_new_down.index[tmp_new_down[s].isin(['0|1', '1|0'])]
                    if len(idx_lst) != 0:
                        idx.append(idx_lst[-1])
                        new = tmp_new_down
                        old = tmp_old_down

                multiplier += 1

            if len(idx) != 0:
                idx = idx[0]
                if new.loc[idx, s][::-1] == old.loc[idx, s]:
                    tmp = beagle_hla.loc[0, f'HLA-{gene} 1']
                    beagle_hla.loc[0, f'HLA-{gene} 1'] = beagle_hla.loc[0, f'HLA-{gene} 2']
                    beagle_hla.loc[0, f'HLA-{gene} 2'] = tmp
            else:
                beagle_hla.loc[0, f'HLA-{gene} 1'] = 'N/A'
                beagle_hla.loc[0, f'HLA-{gene} 2'] = 'N/A'

            beagle_hla.to_csv(output.beagle_phased, header = True, index = False, sep = '\t')

samples_oneKG = read_tsv_as_lst("/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/samples_to_phase.tsv")

rule calculate_and_aggregate_phasing_concordance:
    input:
        beagle_phased = expand("results/phasing/HLA_1KG_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{sample}.1KG.tsv", gene = HLA_GENES, sample = samples_oneKG),
        phase_df = expand("results/phasing/phased_dfs/oneKG_{vcf_version}-{filter}-{gene}.tsv", gene = HLA_GENES, allow_missing = True)
    output:
        concordance_df = "results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}.tsv"
    resources: mem = '30G'
    threads: 4
    run: 
        for i, gene in enumerate(HLA_GENES):
            beagle_hla = pd.concat([pd.read_csv(f"results/phasing/HLA_1KG_BEAGLE/tmp/beagle_phased_per_sample/{gene}.{s}.1KG.tsv", sep = '\t') for s in samples_oneKG])
            beagle_hla = beagle_hla.reset_index(drop = True)
            our_hla = pd.read_csv(f"results/phasing/phased_dfs/oneKG_{wildcards.vcf_version}-{wildcards.filter}-{gene}.tsv", sep = '\t')
            
            tmp = calculate_phasing_concordance(beagle_hla, our_hla, gene)
            if i == 0:
                result = tmp
            else:
                result = pd.concat([result, tmp])

        result.to_csv(output.concordance_df, sep = '\t', header = True, index = False)
