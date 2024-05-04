include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

samples_hc = read_tsv_as_lst(config['samples_hc'])
samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_chip = read_tsv_as_lst(config['samples_chip'])
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
chromosome = [i for i in range(1,23)]

PANEL_NAME = config["PANEL_NAME"]
imp_dir = config["imputation_dir"]

rule retain_chip_sites:
    input:
        lc_vcf = f"results/imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/lc.chr{{chr}}.vcf.gz",
        site = temp(f"results/wip_vcfs/{PANEL_NAME}/chip_sites/chr{{chr}}.tsv")
    resources:
        mem = '30G'
    threads: 4
    shell: """
        mkdir -p results/wip_vcfs/{PANEL_NAME}/chip_sites/

        zgrep -v '#' {input.chip_vcf} | cut -f1,2 > {output.site}
        bcftools view -R {output.site} -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """

rule concat_chip_sites_vcfs:
    input:
        lc_vcf = expand(f"results/wip_vcfs/{PANEL_NAME}/chip_sites/lc.chr{{chr}}.vcf.gz", chr = chromosome)
    output:
        concat = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/lc.vcf.gz",
        bgen = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/lc.bgen"
    resources:
        mem = '100G'
    threads: 16
    params:
        qctool = tools['qctool']
    shell: """
        bcftools concat --threads 8 {input.lc_vcf} | bcftools sort -Oz -o {output.concat}

        gunzip -c {output.concat} | grep '#' > {output.concat}.temp1.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t[%GT\t]\n' \
        {output.concat} >> {output.concat}.temp1.vcf
        bcftools sort -Oz -o {output.concat} {output.concat}.temp1.vcf
        tabix {output.concat}
        rm {output.concat}.temp*

        {params.qctool} -g {output.concat} -og {output.bgen} -bgen-bits 8 -bgen-compression zstd 
    """

rule compute_chip_stats:
    input:
	    bgen = rules.concat_chip_sites_vcfs.output.bgen
    output:
        sqlite = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/lc.sqlite"
    params:
        qctool = tools['qctool']
    shell: """
        {params.qctool} \
        -analysis-name "qc:autosomes" \
        -g {input.bgen} \
        -excl-range X:0- -excl-range Y:0- \
        -snp-stats \
        -osnp sqlite://{output.sqlite}:autosomes \
        -sample-stats \
        -osample sqlite://{output.sqlite}:sample_stats

        {params.qctool} \
        -analysis-name "qc:sex_chromosomes" \
        -g {input.bgen} \
        -incl-range X:0- -incl-range Y:0- \
        -snp-stats \
        -osnp sqlite://{output.sqlite}:sex_chromosomes \
    """

rule thin_stats:
    input:
        db = rules.compute_chip_stats.output.sqlite
    output:
        thinned_ok = touch( f"results/wip_vcfs/{PANEL_NAME}/chip_sites/thinned_ok.ok" ),
        tsv = temp( f"results/wip_vcfs/{PANEL_NAME}/chip_sites/included_variants_included.gen" ),
        tmp = temp( f"results/wip_vcfs/{PANEL_NAME}/chip_sites/tmp.gen" )
    params:
        MAC = 5,
        missing = 10,
        inthinnerator = tools['inthinnerator']
    resources:
        mem = '10G'
    shell: r"""
        mkdir -p results/wip_vcfs/oneKG/chip_sites/

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

rule calculate_chip_PC:
    input:
        sqlite = rules.compute_chip_stats.output.sqlite,
        thin = rules.thin_stats.output.thinned_ok,
        bgen = rules.concat_chip_sites_vcfs.output.bgen
    output:
        variants = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/pc_variants_thin_1bp.txt",
        kinship1 = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/chip_kinship_thin_1bp.all.tsv.gz",
        UDUT1 = f"results/wip_vcfs/{PANEL_NAME}/chip_sites/chip_UDUT_thin_1bp.all.tsv.gz"
    params:
        PCs = 20,
        qctool = tools['qctool']
    resources:
        mem = '10G'
    shell: """
        sqlite3 {input.sqlite} \
        "SELECT chromosome FROM thin_1bp View WHERE result == 'picked'" > {output.variants}
        
        {params.qctool} \
        -analysis-name "PCs:thin_1bp:all" \
        -g {input.bgen} \
        -incl-rsids {output.variants} \
        -kinship {output.kinship1} \
        -UDUT {output.UDUT1} \
        -PCs {params.PCs} \
        -osample sqlite://{input.sqlite}:PCs
    """

rule filter_lc_info:
    input:
        lc_vcf = f"results/imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        info = config['info_filter']
    shell: """
        mkdir -p results/wip_vcfs/{PANEL_NAME}/high_info/

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """

rule filter_lc_maf:
    input:
        lc_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info/lc.chr{{chr}}.vcf.gz",
        af = "data/oneKG_MAFs/oneKG_AF_AFR_chr{chr}.txt"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        info = config['info_filter'],
        maf = config['maf_filter'],
        panel = PANEL_NAME,
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        imp_vcf = input.lc_vcf
        af_txt = input.af

        lc = lcwgsus.read_vcf(imp_vcf).sort_values(by=['chr', 'pos'])
        metadata = lcwgsus.read_metadata(imp_vcf)
        af = lcwgsus.read_af(af_txt)

        lc_af = pd.merge(lc, af, on = common_cols)
        lc_af = lc_af[lc_af['MAF'] > float(params.maf)]
        lc_af = lc_af.drop(columns = 'MAF')

        lcwgsus.save_vcf(lc_af,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )

# The omni5m manifest has col3,4 to be chr and pos
rule prepare_chip_manifest:
    input:
        bpm = config['dense_bpm'],
        egt = config['dense_egt']
    output:
        csv = "data/chip/omni5m/omni5m.csv",
        pos = "data/chip/omni5m/omni5m_sites.tsv"
    resources:
        mem = '30G'
    params:
        picard = tools['picard']
    shell: """
        mkdir -p results/data/chip/omni5m/

        {params.picard} BpmToNormalizationManifestCsv \
        -I {input.bpm} \
        -CF {input.egt} \
        -O {output.csv}

        cut -d ',' -f3,4 {output.csv} | \
        sed 's/,/\t/g' | \
        tail -n +2 > {output.pos}
    """

rule filter_lc_sites:
    input:
        vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af/lc.chr{{chr}}.vcf.gz",
        sites = rules.prepare_chip_manifest.output.pos
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/high_info_high_af_chip_sites/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        panel = PANEL_NAME,
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        imp_vcf = input.vcf
        chip_sites = input.sites

        lc = lcwgsus.read_vcf(imp_vcf).sort_values(by=['chr', 'pos'])
        metadata = lcwgsus.read_metadata(imp_vcf)

        lc = lc.apply(lcwgsus.convert_to_chip_format, axis = 1)
        
        sites = pd.read_table(chip_sites, sep = '\t', names = common_cols, dtype = {'chr': str, 'pos': int}).drop_duplicates(ignore_index = True)
        sites = sites[sites['chr'] == str(wildcards.chr)]
        sites['chr'] = sites['chr'].astype(int)

        lc_sites = pd.merge(lc, sites, on = common_cols)

        lcwgsus.save_vcf(lc_sites,
             metadata,
             prefix='chr',
             outdir="results/wip_vcfs/" + params.panel + "/high_info_high_af_chip_sites/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )