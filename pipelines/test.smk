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
        # v = "/well/band/users/rbx225/test_files/GAM013489_oneKG/chr3.vcf.gz",
        t = "test.txt"

rule test:
    output:
        # v = "/well/band/users/rbx225/test_files/GAM013489_oneKG/chr3.vcf.gz",
        t = "test.txt"
    params:
        imputation_dir = config['imputation_dir']
    run:
        x = params.imputation_dir.split('/')[-2]
        print(x)

        shell("echo test.txt")