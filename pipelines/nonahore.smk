configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append('/well/band/users/rbx225/software/QUILT_test/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *

sys.path.append('/well/band/users/rbx225/GAMCC/scripts/lcSV/')
sys.path.append('/Users/sus_zhang/Desktop/Suuuuuuuus/Low Coverage Data/gamcc/scripts/lcSV/')
from lcSV import *

replicates = 100

rule all:
    input:
        pickle = expand('results/nonahore/simulate/plausibility/rep{rep}/eval.pickle', rep = [i for i in range(replicates)])

rule simulate_nonahore:
    output:
        pickle = 'results/nonahore/simulate/plausibility/rep{rep}/eval.pickle'
    threads: 16
    params:
        odir = 'results/nonahore/simulate/plausibility/rep{rep}/'
    script:
        '/well/band/users/rbx225/GAMCC/scripts/simulate_nonahore.py'