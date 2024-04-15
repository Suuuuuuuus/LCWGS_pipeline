configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

imp_dir = config['imputation_dir']

x = ['a.txt', 'b.txt']
y = ['1.txt', '2.txt']
z = ['s.txt', 't.txt']

analysis = [i for i in range(len(x))]

rule all:
    input:
        res = expand(z[int("{n}")], n = analysis)

rule test:
    input:
        a = x[{n}],
        b = y[{n}]
    output:
        t = z[{n}]
    run:
        print(output.t)