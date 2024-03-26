include: "reference.smk"
include: "chip_qc.smk"
#include: "chip_imputation.smk"

#include: "test.smk"
include: "auxiliary.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_chip = read_tsv_as_lst(config['samples_chip'])
chromosome = [i for i in range(1,23)]

# The followings are global parameters:
clean_fastq = config['clean_fastq']
reheader = config['reheader']
concatenate = config['concatenate']

rule reference_all:
    input:
        amb = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.amb" if concatenate else "data/references/GRCh38.fa.amb",
        ann = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.ann" if concatenate else "data/references/GRCh38.fa.ann",
        bwt = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.bwt" if concatenate else "data/references/GRCh38.fa.bwt",
        pac = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.pac" if concatenate else "data/references/GRCh38.fa.pac",
        sa = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.sa" if concatenate else "data/references/GRCh38.fa.sa"

# chip_thinning = ['thin_1bp', 'thin_50kb', 'thin_100kb']
chip_thinning = ['thin_1bp']

# The following rule has to be run at least twice:
# The first time one needs to generate the qc files. Then one should use the provided jupyter notebook or Python or R script to generate sites/sample information to be retained or removed.
# The second round one could get cleaned vcf files.
rule chip_qc_all:
    input:
        chip_vcf = "results/chip/vcf/chip_genotype.vcf.gz",
        chip_samples = "results/chip/vcf/chip_genotype.sample",
        chip_bgen = "results/chip/bgen/chip.bgen",
        chip_stats = "results/chip/qc/chip.qc.sqlite",
        variants = expand("results/chip/qc/PCs/pc_variants_{thinning}.txt", thinning = chip_thinning),
        kinship1 = expand("results/chip/qc/PCs/chip_kinship_{thinning}.all.tsv.gz", thinning = chip_thinning),
        UDUT1 = expand("results/chip/qc/PCs/chip_UDUT_{thinning}.all.tsv.gz", thinning = chip_thinning),
        #kinship2 = expand("results/chip/qc/PCs/chip_kinship_{thinning}.exclude-duplicates.tsv.gz", thinning = chip_thinning),
        #UDUT2 = expand("results/chip/qc/PCs/chip_UDUT_{thinning}.exclude-duplicates.tsv.gz", thinning = chip_thinning)
        #vcf_qced = "results/chip/vcf/chip_qced.vcf.gz"

rule chip_imputation_all:
    input:

rule test_all:
    input:
        vcf = "results/tmp/{id}.{chr}.txt"
