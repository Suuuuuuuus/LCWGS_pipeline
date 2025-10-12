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
from hla_phase_functions import *
from hla_align_functions import *

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
SEED = 42

rule all:
    input:
        imputed_merged_ref_v3570 = expand("results/hla/imputation/QUILT_HLA_result_db/{id}/{hla_gene}/quilt.hla.output.combined.all.txt", hla_gene = ['A'], id = ['IDT0481'])



rule hla_imputation_db:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/HLA{hla_gene}fullallelesfilledin.RData"
    output:
        bamlist = temp("results/hla/imputation/QUILT_HLA_result_db/{id}/{id}.{hla_gene}.tsv"),
        imputed = "results/hla/imputation/QUILT_HLA_result_db/{id}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '20G'
    threads: 2
    conda: "sus"
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_db/",
        seed = SEED
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_optimal/{wildcards.id}/{wildcards.hla_gene}/
        echo {input.bam} >> {output.bamlist}

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_db/{wildcards.id}/{wildcards.hla_gene}/" \
        --quilt_seed={params.seed} \
        --bamlist={output.bamlist} \
        --region={wildcards.hla_gene} \
        --n_seek_iterations=2 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """