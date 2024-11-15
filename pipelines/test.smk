configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

rule all:
    input:
        res1 = "results/hla/imputation/QUILT_HLA_result_method/IDT04811/DRB1/quilt.hla.output.combined.all.txt",
        res2 = "results/hla/imputation/QUILT_HLA_result_method/IDT04812/DRB1/quilt.hla.output.combined.all.txt",
        res3 = "results/hla/imputation/QUILT_HLA_result_method/IDT04813/DRB1/quilt.hla.output.combined.all.txt",
        res4 = "results/hla/imputation/QUILT_HLA_result_method/IDT04814/DRB1/quilt.hla.output.combined.all.txt"
        
rule imp1:
    input:
        bam = "data/bams/IDT0481.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLADRB1fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/DRB1.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/IDT04811.DRB1.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/IDT04811/DRB1/quilt.hla.output.combined.all.txt"
    resources:
        mem = '80G'
    threads: 6
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/IDT04811/DRB1/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/IDT04811/DRB1/" \
        --bamlist={output.bamlist} \
        --region=DRB1 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.DRB1.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

rule imp2:
    input:
        bam = "data/bams/IDT0481.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLADRB1fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/DRB1.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/IDT04812.DRB1.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/IDT04812/DRB1/quilt.hla.output.combined.all.txt"
    resources:
        mem = '120G'
    threads: 16
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/IDT04812/DRB1/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/IDT04812/DRB1/" \
        --bamlist={output.bamlist} \
        --region=DRB1 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.DRB1.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

rule imp3:
    input:
        bam = "data/bams/IDT0481.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLADRB1fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/DRB1.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/IDT04813.DRB1.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/IDT04813/DRB1/quilt.hla.output.combined.all.txt"
    resources:
        mem = '80G'
    threads: 6
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus2"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/IDT04813/DRB1/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/IDT04813/DRB1/" \
        --bamlist={output.bamlist} \
        --region=DRB1 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.DRB1.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

rule imp4:
    input:
        bam = "data/bams/IDT0481.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLADRB1fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/DRB1.ssv'
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/IDT04814.DRB1.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/IDT04814/DRB1/quilt.hla.output.combined.all.txt"
    resources:
        mem = '120G'
    threads: 16
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus2"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/IDT04814/DRB1/
        ls {input.bam} > {output.bamlist}
        ulimit -n 50000

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/IDT04814/DRB1/" \
        --bamlist={output.bamlist} \
        --region=DRB1 \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.DRB1.haplotypes.RData \
        --dict_file={params.fa_dict}
    """