configfile: "pipelines/config.json"
include: "auxiliary.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import lcwgSus

chromosome = [i for i in range(1,23)]

samples_hc = list(pd.read_table(config['samples_hc'], header = None, names = ['Code'])['Code'].values)
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
test_hc = ids_1x_all[:2]

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
        script = "scripts/convert_chip.R"
    shell: """
        mkdir -p results/chip/vcf/
        
        Rscript --vanilla {params.script} \
        --manifest {input.manifest} \
        --genotypes {input.genotypes} \
        --samples {input.samples} \
        --output {output.tmp}

        bgzip {output.tmp}
        
        mkdir -p results/chip/bgen/
        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
        -g {output.vcf} -s {output.samples} -og {output.bgen} -bgen-bits 8 -bgen-compression zstd
    """	

# Compute SNP and sample stats across autosomes and sex chromosomes separately
rule compute_chip_stats:
    input:
	    bgen = rules.genotype_chip.output.bgen,
	    samples = rules.genotype_chip.output.samples
    output:
        sqlite = "results/chip/qc/chip.qc.sqlite"
    shell: """
        mkdir -p results/chip/qc/

        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
        -analysis-name "qc:autosomes" \
        -g {input.bgen} \
        -s {input.samples} \
        -excl-range X:0- -excl-range Y:0- \
        -snp-stats \
        -osnp sqlite://{output.sqlite}:autosomes \
        -sample-stats \
        -osample sqlite://{output.sqlite}:sample_stats

        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
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
        missing = 10
    resources:
        mem = '10G'
    shell: r"""
        mkdir -p results/chip/qc/PCs/

        sqlite3 -header -separator $'\t' {input.db} \
        "SELECT rsid AS SNPID, rsid, chromosome, position, alleleA, alleleB FROM autosomesView WHERE (alleleA_count >= {params.MAC}) AND (alleleB_count >= {params.MAC}) AND \`NULL\` < {params.missing}" > {output.tmp}

        tail -n +2 {output.tmp} > {output.tsv}

        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/inthinnerator_v2.2.2 \
        -analysis-name thin_100kb \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 100kb \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_100kb

        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/inthinnerator_v2.2.2 \
        -analysis-name thin_5kb \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 50kb \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_50kb
  
        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/inthinnerator_v2.2.2 \
        -analysis-name thin_1bp \
        -g {output.tsv} \
        -suppress-excluded \
        -min-distance 1bp \
        -excl-range 06:25000000-40000000 \
        -o sqlite://{input.db}:thin_1bp
    """

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
        PCs = 20
    resources:
        mem = '10G'
    shell: """
        sqlite3 {input.sqlite} \
        "SELECT chromosome FROM {wildcards.thinning}View WHERE result == 'picked'" > {output.variants}
        
        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
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
        /well/band/users/rbx225/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
        -analysis-name "PCs:{wildcards.thinning}:exclude-duplicates" \
        -g {input.bgen} -s {input.samples} \
        -incl-rsids {input.variants} \
        -kinship {output.kinship2} \
        -UDUT {output.UDUT2} \
        -PCs {params.PCs} \
        -excl-samples-where "ID = 'GAM370894'" \
        -excl-samples-where "ID = 'GAM916387'" \
        -excl-samples-where "ID = 'GAM654203'" \
        -osample sqlite://{input.sqlite}:PCs
    """
    