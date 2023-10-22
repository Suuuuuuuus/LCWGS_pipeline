configfile: "pipelines/config.json"

import pandas as pd
import sys
sys.path.append("scripts")
import lcwgSus

config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
chromosome = [i for i in range(1,23)]

rule test:
    input:
        ss_bam = expand("data/subsampled_bams/{id}_subsampled.bam", id = ids_1x_all)
    output:
        graph = "results/test.txt"
    params:
        dedup = config["dedup"]
    run:
        with open(output.graph, "w") as f:
            sys.stdout = f
            print(type(input.ss_bam))
            print()
            for i in input.ss_bam:
                print(i)
