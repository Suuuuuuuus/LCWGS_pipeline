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
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

rule beagle_phase_1kg:
    input:
        vcf = "results/phasing/HLA_1KG_BEAGLE/unphased.1KG.chr6.vcf.gz"
    output:
        phased_vcf = "results/phasing/HLA_1KG_BEAGLE/phased.1KG.chr6.vcf.gz"
    params:
        beagle = tools['beagle'],
        output_prefix = "results/phasing/HLA_1KG_BEAGLE/phased.1KG.chr6"
    shell: """
        mkdir -p results/phasing/HLA_1KG_BEAGLE/
        {params.beagle} gt={input.vcf} out={params.output_prefix}
    """
