configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

samples_hc = read_tsv_as_lst(config['samples_hc'])

# This script is borrowed from Dr Gavin Band. The original script is capable of coping with multiple builds (GRCh37/38) by providing alternative manifest files. In this pipeline, the behavior is temporarily disabled.

'''
`.annot.db` file contains annotation information sent from the company.
`Genotype` file is the actual chip result in raw format.
`Sample_Table` file contains sample information.
'''

# Using annotation file to convert chip genotypes from raw format to vcf file
rule genotype_chip:
    input:
        genotypes = config["chip_genotypes"],
        samples = config["chip_samples"],
        manifest = config["chip_annotation"]
    output:
        tmp = temp("results/chip/vcf/chip_genotype.vcf"),
        vcf = "results/chip/vcf/chip_genotype.vcf.gz",
        samples = "results/chip/vcf/chip_genotype.sample",
        bgen = "results/chip/bgen/chip.bgen"
    params:
        script = "scripts/convert_chip.R",
        qctool = tools['qctool']
    shell: """
        mkdir -p results/chip/vcf/
        
        Rscript --vanilla {params.script} \
        --manifest {input.manifest} \
        --genotypes {input.genotypes} \
        --samples {input.samples} \
        --output {output.tmp}

        bgzip {output.tmp}
        
        mkdir -p results/chip/bgen/
        {params.qctool} \
        -g {output.vcf} -s {output.samples} -og {output.bgen} -bgen-bits 8 -bgen-compression zstd
    """	

# Compute SNP and sample stats across autosomes and sex chromosomes separately
rule compute_chip_stats:
    input:
	    bgen = rules.genotype_chip.output.bgen,
	    samples = rules.genotype_chip.output.samples
    output:
        sqlite = "results/chip/qc/chip.qc.sqlite"
    params:
        qctool = tools['qctool']
    shell: """
        mkdir -p results/chip/qc/

        {params.qctool} \
        -analysis-name "qc:autosomes" \
        -g {input.bgen} \
        -s {input.samples} \
        -excl-range X:0- -excl-range Y:0- \
        -snp-stats \
        -osnp sqlite://{output.sqlite}:autosomes \
        -sample-stats \
        -osample sqlite://{output.sqlite}:sample_stats

        {params.qctool} \
        -analysis-name "qc:sex_chromosomes" \
        -g {input.bgen} \
        -s {input.samples} \
        -incl-range X:0- -incl-range Y:0- \
        -snp-stats \
        -osnp sqlite://{output.sqlite}:sex_chromosomes \
    """

rule thin_stats:
    input:
        db = rules.compute_chip_stats.output.sqlite
    output:
        thinned_ok = touch( "results/chip/qc/PCs/thinned_ok.ok" ),
        tsv = temp( "results/chip/qc/PCs/included_variants_included.gen" ),
        tmp = temp( "results/chip/qc/PCs/tmp.gen" )
    params:
        MAC = 5,
        missing = 10,
        inthinnerator = tools['inthinnerator']
    resources:
        mem = '10G'
    shell: r"""
        mkdir -p results/chip/qc/PCs/

        sqlite3 -header -separator $'\t' {input.db} \
        "SELECT rsid AS SNPID, rsid, chromosome, position, alleleA, alleleB FROM autosomesView WHERE (alleleA_count >= {params.MAC}) AND (alleleB_count >= {params.MAC}) AND \`NULL\` < {params.missing}" > {output.tmp}

        tail -n +2 {output.tmp} > {output.tsv}

        {params.inthinnerator} \
        -analysis-name thin_1bp \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 1bp \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_1bp
    """

'''
        {params.inthinnerator} \
        -analysis-name thin_100kb \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 100kb \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_100kb

        {params.inthinnerator} \
        -analysis-name thin_5kb \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 50kb \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_50kb
  
'''

rule calculate_chip_PC:
    input:
        sqlite = rules.compute_chip_stats.output.sqlite,
        thin = rules.thin_stats.output.thinned_ok,
        bgen = rules.genotype_chip.output.bgen,
        samples = rules.genotype_chip.output.samples
    output:
        variants = "results/chip/qc/PCs/pc_variants_{thinning}.txt",
        kinship1 = "results/chip/qc/PCs/chip_kinship_{thinning}.all.tsv.gz",
        UDUT1 = "results/chip/qc/PCs/chip_UDUT_{thinning}.all.tsv.gz"
    params:
        PCs = 20,
        qctool = tools['qctool']
    resources:
        mem = '10G'
    shell: """
        sqlite3 {input.sqlite} \
        "SELECT chromosome FROM {wildcards.thinning}View WHERE result == 'picked'" > {output.variants}
        
        {params.qctool} \
        -analysis-name "PCs:{wildcards.thinning}:all" \
        -g {input.bgen} -s {input.samples} \
        -incl-rsids {output.variants} \
        -kinship {output.kinship1} \
        -UDUT {output.UDUT1} \
        -PCs {params.PCs} \
        -osample sqlite://{input.sqlite}:PCs
    """
# The previous rule should be `SELECT rsid` but somehow the order of columns are messed up...

# This rule has to be run after the previous and need some manual input from the qc metrics. This behavior should be adjusted later.
rule exclude_chip_dup_samples:
    input:
        sqlite = rules.compute_chip_stats.output.sqlite,
        bgen = rules.genotype_chip.output.bgen,
        samples = rules.genotype_chip.output.samples,
        variants = rules.calculate_chip_PC.output.variants
    output:
        kinship2 = "results/chip/qc/PCs/chip_kinship_{thinning}.exclude-duplicates.tsv.gz",
        UDUT2 = "results/chip/qc/PCs/chip_UDUT_{thinning}.exclude-duplicates.tsv.gz"
    params:
        PCs = 20
    resources:
        mem = '10G'
    shell: """
        {params.qctool} \
        -analysis-name "PCs:{wildcards.thinning}:exclude-duplicates" \
        -g {input.bgen} -s {input.samples} \
        -incl-rsids {input.variants} \
        -kinship {output.kinship2} \
        -UDUT {output.UDUT2} \
        -PCs {params.PCs} \
        -excl-samples-where "ID = 'GAM370894'" \
        -excl-samples-where "ID = 'GAM916387'" \
        -osample sqlite://{input.sqlite}:PCs
    """
   
# Try if -T and -R make any difference (it shouldn't)
# This rule has to be run separately with the previous, as it needs manual creation of the drop_samples and retain_sites file from the qc results. See the provided jupyter notebook.
rule clean_chip_vcf:
    input:
        vcf = rules.genotype_chip.output.vcf
    output:
        vcf_qced = "results/chip/vcf/chip_qced.vcf.gz"
    params:
        retain_sites = "results/chip/vcf/retain_sites.tsv"
    resources:
        mem = '20G'
    shell: """
        tabix {input.vcf}

        if [[ -d {params.retain_sites}]]
        then
            bcftools view -T {params.retain_sites} -Oz -o {output.vcf_qced} {input.vcf}
            tabix {output.vcf_qced}
        else
            touch {output.vcf_qced} # A place holder to avoid errors
        fi
    """

rule calculate_PCA:
    input:
        vcf = "results/chip/vcf/chip_qced.vcf.gz"
    output:
        bed = temp("results/chip/qc/PCs/chip_pca.bed"),
        bim = temp("results/chip/qc/PCs/chip_pca.bim"),
        fam = temp("results/chip/qc/PCs/chip_pca.fam"),
        PC = "results/chip/qc/PCs/PCs.eigenvec",
        tmp_vcf = temp("results/chip/qc/PCs/tmp.vcf.gz"),
        name_txt = temp("results/chip/qc/PCs/name.txt")
    params:
        num_PCs = 10,
        plink_name = "results/chip/qc/PCs/chip_pca",
        PC_name = "results/chip/qc/PCs/PCs",
        gam_to_exclude = "data/sample_tsvs/PC_exclude/gam_to_exclude.tsv"
    resources:
        mem = '10G'
    shell: """
        cp {params.gam_to_exclude} {output.name_txt}
        bcftools view -S ^{output.name_txt} -Oz -o {output.tmp_vcf} {input.vcf}
        plink --vcf {output.tmp_vcf} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """