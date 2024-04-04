configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

rule all:
    input:
        "test.txt"

rule test:
    input:
        "results/imputation/vcfs/oneKG/"
    output:
        "test.txt"
    shell: """

        for i in {{1..10}}; do echo $i >> {output}; done
    """