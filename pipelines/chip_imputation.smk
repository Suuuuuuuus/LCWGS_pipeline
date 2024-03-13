configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
test_hc = samples_lc[:2]

hc_panel = config["hc_panel"]
