configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

imp_dir = config['imputation_dir']

rule all:
    input:
        v = "/well/band/users/rbx225/test_files/GAM013489_oneKG/chr1.vcf.gz"

rule test:
    output:
        v = "/well/band/users/rbx225/test_files/GAM013489_oneKG/chr1.vcf.gz"
    run:
        print('hi')
        
        shell("""
            gunzip {vcf}; bgzip /well/band/users/rbx225/{panel}/GAM013489_oneKG/chr{c}.vcf
        """.format(vcf = input.v, panel = "test_files", c = "1"))