include: "imputation_calculation.smk"

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

samples_lc = read_tsv_as_lst(config['samples_lc'])
panels = config['panels']

chromosome = [i for i in range(1,23)]

rule lc_chip_all:
    input:
        h_report = expand("results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.tsv", chr = chromosome),
        h_impacc = expand("results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv", chr = chromosome),
        v_report = expand("results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.tsv", chr = chromosome),
        v_impacc = expand("results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv", chr = chromosome),
        r2NRC_h = "graphs/imputation/lc_chip/all_samples/by_variant/r2_NRC.png",
        ccd_h = "graphs/imputation/lc_chip/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = "graphs/imputation/lc_chip/all_samples/by_sample/r2_NRC.png",
        ccd_v = "graphs/imputation/lc_chip/all_samples/by_sample/ccd_by_genotype.png"

rule manipulate_vcf:
    input:
        info