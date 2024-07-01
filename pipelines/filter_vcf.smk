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
from lcwgsus.variables import *

chromosome = [i for i in range(1,23)]

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

rule retain_chip_sites:
    input:
        lc_vcf = "results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz"
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/chip_sites/lc.chr{chr}.vcf.gz",
        site = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/chr{chr}.tsv")
    resources:
        mem = '30G'
    threads: 4
    shell: """
        mkdir -p results/wip_vcfs/{wildcards.panel}/vanilla/chip_sites/

        zgrep -v '#' {input.chip_vcf} | cut -f1,2 > {output.site}
        bcftools view -R {output.site} -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """
'''
rule concat_chip_sites_vcfs:
    input:
        lc_vcf = expand("results/wip_vcfs/{panel}/vanilla/chip_sites/lc.chr{chr}.vcf.gz", chr = chromosome)
    output:
        concat = "results/wip_vcfs/{panel}/vanilla/chip_sites/lc.vcf.gz"
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
    """

rule calculate_PCA:
    input:
        vcf = rules.concat_chip_sites_vcfs.output.concat
    output:
        bed = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/lc_pca.bed"),
        bim = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/lc_pca.bim"),
        fam = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/lc_pca.fam"),
        PC = "results/wip_vcfs/{panel}/vanilla/chip_sites/PCs.eigenvec",
        name_txt = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/name.txt"),
        tmp_vcf = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/tmp.vcf.gz"),
        tmp_vcf1 = temp("results/wip_vcfs/{panel}/vanilla/chip_sites/tmp1.vcf.gz")
    params:
        num_PCs = 10,
        plink_name = "results/wip_vcfs/{panel}/vanilla/chip_sites/lc_pca",
        PC_name = "results/wip_vcfs/{panel}/vanilla/chip_sites/PCs",
        chip_gm_names = "data/sample_tsvs/chip_gm_names.tsv"
    resources:
        mem = '10G'
    shell: """
        cp {params.chip_gm_names} {output.name_txt}
        bcftools view -S {output.name_txt} -Oz -o {output.tmp_vcf} {input.vcf}
        cp {params.gm_to_exclude} {output.name_txt}
        bcftools view -S ^{output.name_txt} -Oz -o {output.tmp_vcf1} {output.tmp_vcf}
        plink --vcf {output.tmp_vcf1} --make-bed --out {params.plink_name}
        plink --bfile {params.plink_name} --pca {params.num_PCs} --out {params.PC_name}
    """
'''
rule filter_lc_info:
    input:
        lc_vcf = "results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz"
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info/lc.chr{chr}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        info = config['info_filter']
    shell: """
        mkdir -p results/wip_vcfs/{wildcards.panel}/vanilla/high_info/

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """

rule filter_lc_maf:
    input:
        lc_vcf = "results/wip_vcfs/{panel}/vanilla/high_info/lc.chr{chr}.vcf.gz",
        af = "data/oneKG_MAFs/oneKG_AF_AFR_chr{chr}.txt"
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        info = config['info_filter'],
        maf = config['maf_filter'],
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos', 'ref', 'alt']

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
            outdir="results/wip_vcfs/" + wildcards.panel + "/vanilla/high_info_high_af/",
            save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
            )

rule filter_lc_sites:
    input:
        vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz",
        sites = rules.prepare_chip_manifest.output.pos
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_chip_sites/lc.chr{chr}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos']

        imp_vcf = input.vcf
        chip_sites = input.sites

        lc = lcwgsus.read_vcf(imp_vcf).sort_values(by = ['chr', 'pos'])
        metadata = lcwgsus.read_metadata(imp_vcf)

        lc = lc.apply(lcwgsus.convert_to_chip_format, axis = 1)
        
        sites = pd.read_table(chip_sites, sep = '\t', names = common_cols, dtype = {'chr': str, 'pos': int}).drop_duplicates(ignore_index = True)
        sites = sites[sites['chr'] == str(wildcards.chr)]
        sites['chr'] = sites['chr'].astype(int)

        lc_sites = pd.merge(lc, sites, on = common_cols)

        lcwgsus.save_vcf(lc_sites,
            metadata,
            prefix='chr',
            outdir="results/wip_vcfs/" + wildcards.panel + "/vanilla/high_info_high_af_chip_sites/",
            save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
            )

rule filter_oneKG_low_confidence_regions:
    input:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz"
    output:
        filtered_vcf = temp("results/wip_vcfs/{panel}/vanilla/high_info_high_af_high_conf/lc.chr{chr}.vcf"),
        filtered_vcf_gz = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_high_conf/lc.chr{chr}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        bed = "data/bedgraph/strict.all.bed"
    shell: """
        mkdir -p results/wip_vcfs/{wildcards.panel}/vanilla/high_info_high_af_high_conf/

        gunzip -c {input.filtered_vcf} | grep '#' > {output.filtered_vcf}
        tabix -R {params.bed} {input.filtered_vcf} >> {output.filtered_vcf}

        bgzip -f {output.filtered_vcf}
        tabix {output.filtered_vcf_gz}
        touch {output.filtered_vcf}
    """
'''
rule filter_giab_low_confidence_regions:
    input:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz"
    output:
        filtered_vcf = temp("results/wip_vcfs/{panel}/vanilla/high_info_high_af_giab_high_conf/lc.chr{chr}.vcf"),
        filtered_vcf_gz = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_giab_high_conf/lc.chr{chr}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        bed = "data/bedgraph/GIAB_not_in_difficult_regions.bed"
    shell: """
        mkdir -p results/wip_vcfs/{panel}/vanilla/high_info_high_af_giab_high_conf/

        gunzip -c {input.filtered_vcf} | grep '#' > {output.filtered_vcf}
        tabix -R {params.bed} {input.filtered_vcf} >> {output.filtered_vcf}
        bgzip {output.filtered_vcf}
        tabix {output.filtered_vcf_gz}
        touch {output.filtered_vcf}
    """

rule filter_giab_lc_sites:
    input:
        vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_giab_high_conf/lc.chr{chr}.vcf.gz",
        sites = rules.prepare_chip_manifest.output.pos
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_giab_high_conf_chip_sites/lc.chr{chr}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        panel = panel,
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
             outdir="results/wip_vcfs/" + params.panel + "/vanilla/high_info_high_af_giab_high_conf_chip_sites/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )
'''
rule filter_oneKG_lc_sites:
    input:
        vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_high_conf/lc.chr{chr}.vcf.gz",
        sites = rules.prepare_chip_manifest.output.pos
    output:
        filtered_vcf = "results/wip_vcfs/{panel}/vanilla/high_info_high_af_high_conf_chip_sites/lc.chr{chr}.vcf.gz"
    resources:
        mem = '60G'
    threads: 8
    params:
        chrom = "{chr}"
    run:
        common_cols = ['chr', 'pos']

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
            outdir="results/wip_vcfs/" + wildcards.panel + "/vanilla/high_info_high_af_high_conf_chip_sites/",
            save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
            )