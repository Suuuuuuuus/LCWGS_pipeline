configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
hc_panel = config["hc_panel"]
chromosome = [i for i in range(1,23)]
variant_types = ['snps', 'indels']
concatenate = config["concatenate"]

