import io
import os
import re
import sys
import csv
import gzip
import subprocess
import pandas as pd
import numpy as np
import itertools
import collections
sys.path.append('/well/band/users/rbx225/software/lcwgsus/')
import lcwgsus
from lcwgsus.variables import *

hla_gene_information = pd.read_csv('/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv', sep = ' ')
samples = lcwgsus.read_tsv_as_lst('/well/band/users/rbx225/GAMCC/data/sample_tsvs/fv_gm_names.tsv')

hla = lcwgsus.read_hla_direct_sequencing(retain = 'fv')
hla = hla.drop(columns = ['One field1', 'One field2'])

def keep_first_allele(a):
    if '/' in a:
        alleles = a.split('/')
        onefield = alleles[0].split(':')[0]
        second_fields = [
            int(re.match(r"(\d+)", allele.split(':')[1]).group(1)) 
            for allele in alleles
        ]
        smallest_second_field = min(second_fields)
        return f"{onefield}:{smallest_second_field:02d}"
    elif 'XX' in a:
        return '-9'
    else:
        return a

hla['Two field1'] = hla['Two field1'].apply(keep_first_allele)
hla['Two field2'] = hla['Two field2'].apply(keep_first_allele)

colnames = ['Sample ID'] + [label for g in HLA_GENES for label in [f'HLA-{g} 1', f'HLA-{g} 2']]
hlatypes = pd.DataFrame(columns = colnames)
for s in hla['SampleID'].unique():
    tmp = hla[hla['SampleID'] == s]
    row = [s] + tmp[['Two field1', 'Two field2']].values.ravel().tolist()
    hlatypes.loc[len(hlatypes)] = row
    
sl = pd.read_csv(SAMPLE_LINKER_FILE)
sl = sl[~sl['Sample_Name'].str.contains('mini')]
sl = {k:v for k, v in zip(sl['Chip_Name'], sl['Sample_Name'])}
hlatypes['Sample ID'] = hlatypes['Sample ID'].apply(lambda x: sl[x])

phased_vcf = "/well/band/users/rbx225/GAMCC/results/imputation/vcfs/malariaGen_v1_b38/quilt.chr6.vcf.gz"
# phased_vcf = "/well/band/users/rbx225/GAMCC/results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/emp_vcf/chr6.empiricalDose.vcf.gz"

start = 25000000
end = 34000000

command = f"bcftools view -s {','.join(samples)} -r chr6:" + str(int(start)) + "-" + str(int(end)) + " " + phased_vcf
vcf = subprocess.run(command, shell = True, capture_output = True, text = True).stdout[:-1].split('\n')
vcf = [i.split('\t') for i in vcf if '##' not in i]
vcf = pd.DataFrame(vcf)
vcf.columns = vcf.iloc[0]  # Set the first row as the header
vcf = vcf[1:].reset_index(drop = True)
vcf['POS'] = vcf['POS'].astype(int)
vcf = vcf[vcf.columns[:9].tolist() + hlatypes['Sample ID'].tolist()]
# vcf = vcf[vcf.columns[:9] + hlatypes['Sample ID'].tolist()]

distinct_alleles = {}
for g in HLA_GENES:
    alleles = np.append(hlatypes[f'HLA-{g} 1'].astype(str).to_numpy(), hlatypes[f'HLA-{g} 2'].astype(str).to_numpy())
    alleles = np.unique(alleles)
    alleles = alleles[alleles != 'nan']
    alleles = alleles[alleles != '-9']
    distinct_alleles[g] = alleles

for g in HLA_GENES:
    r = hla_gene_information[hla_gene_information['Name'] == f'HLA-{g}']
    start = r['Start'].values[0]
    all_snp_pos = vcf['POS'].tolist()
    position = start
    n_alleles = 0

    while n_alleles < len(distinct_alleles[g]):
        if position in all_snp_pos:
            pass
        else:
            twofield = distinct_alleles[g][n_alleles]
            common_cols = ['chr6', position, f'HLA_{g}*{twofield}', 'T', 'A', '.', '.', '.', 'GT']
            allele1 = (hlatypes[f'HLA-{g} 1'] == twofield).astype(int).values
            allele2 = (hlatypes[f'HLA-{g} 2'] == twofield).astype(int).values

            gts = [f'{a1}/{a2}' for a1, a2 in zip(allele1, allele2)]
            vcf.loc[len(vcf)] = common_cols + gts
            n_alleles += 1
        position += 1

vcf = vcf.sort_values(by = 'POS', ascending = True).reset_index(drop = True)

metadata = lcwgsus.read_metadata(phased_vcf, filetype = 'gzip', comment = '#', new_cols = None)
metadata[-1] = '\t'.join(vcf.columns) + '\n'
lcwgsus.save_vcf(vcf,
             metadata,
             rezip = True,
             prefix='',
             outdir='/well/band/users/rbx225/GAMCC/results/phasing/HLA_GAMCC_BEAGLE/',
             save_name='unphased.GAMCC.chr6.vcf.gz'
             )