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

chromosome = [i for i in range(1,23)]

QUILT_HOME = config["QUILT_HOME"]
mini_analysis_dir = config["mini_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

samples_lc = read_tsv_as_lst(config['samples_lc'])

rule mini_clean_bam:
    input:
        bam = "data/subsampled_bams/{id}_subsampled.bam"
    output:
        bam = "data/mini_bams/{id}.bam",
        bai = "data/mini_bams/{id}.bam.bai",
        tmp1 = temp("data/mini_bams/{id}.tmp1.bam"),
        metric = temp("data/mini_bams/{id}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/mini_bams/tmp/{id}/",
        sample = "{id}"
    shell: """
        mkdir -p {params.tmpdir}

        samtools index {input.bam}

        picard AddOrReplaceReadGroups \
        -VERBOSITY ERROR \
        -I {input.bam} \
        -O {output.tmp1} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}

        picard FixMateInformation -I {output.tmp1}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

rule prepare_bamlist:
    input:
        bams = expand("data/mini_bams/{id}.bam", id = samples_lc)
    output:
        bamlist = "results/mini_imputation/bamlist.txt"
    shell: """
        mkdir -p {mini_analysis_dir}
        
        ls data/mini_bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {mini_analysis_dir} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/mini_imputation/refs/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME,
        threads=1
    shell: """
        mkdir -p results/mini_imputation/refs/
        bcftools norm -m+ {input.vcf} | bcftools view -Oz -o {output.tmp_vcf} -m2 -M2 -v snps

        tabix {output.tmp_vcf}

        bcftools convert --haplegendsample results/mini_imputation/refs/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/quilt_accessories/determine_chunks.R"
    output:
        json = "results/mini_imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {mini_analysis_dir:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """

rule prepare_ref:
    input:
        json = "results/mini_imputation/regions.json",
        hap = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/mini_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        recomb = f"results/mini_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/mini_imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem_mb = 30000
    params:
        threads = 8
    shell: """
        mkdir -p results/mini_imputation/refs/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/mini_imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        nGen={NGEN}, \
        reference_haplotype_file="{input.hap}" ,\
        reference_legend_file="{input.legend}", \
        genetic_map_file="{input.recomb}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        output_file="{output.RData}")'
    """

rule quilt_imputation:
    input:
        bamlist = "results/mini_imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/mini_imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources:
        mem_mb = 30000
    resources:
        mem_mb = 30000
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    params:
        panel = PANEL_NAME
    shell: """
        ## set a seed here, randomly, so can try to reproduce if it fails
        SEED=`echo $RANDOM`
        mkdir -p "results/mini_imputation/vcfs/{params.panel}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/mini_imputation/refs/RData/other/", \
        chr="chr{wildcards.chr}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        bamlist="{input.bamlist}", \
        prepared_reference_filename="{input.RData}", \
        output_filename="{output.vcf}", \
        seed='${{SEED}}')'
    """

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/mini_imputation/regions.json"
if os.path.exists(file):
    with open(file) as json_file:
        REGIONS = json.load(json_file)

vcfs_to_concat={}
final_vcfs=[]
for chr in chromosome:
    start=REGIONS[str(chr)]["start"]
    end=REGIONS[str(chr)]["end"]
    per_chr_vcfs=[]
    for i in range(0, start.__len__()):
        regionStart=start[i]
        regionEnd=end[i]
        file="results/mini_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/mini_imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule concat_quilt_vcf:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/mini_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
    resources:
        mem_mb = 30000
    params:
        threads = 1,
        input_string=get_input_vcfs_as_string
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        if [ -e {output.vcf}.temp1.vcf.gz ]; then
            rm {output.vcf}.temp1.vcf.gz
        fi
        if [ -e {output.vcf}.temp2.vcf.gz ]; then
            rm {output.vcf}.temp2.vcf.gz
        fi

        bcftools concat \
        --ligate-force \
        --output-type z \
        --output {output.vcf}.temp1.vcf.gz \
        {params.input_string}

        gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS\t[%GT:%GP:%DS\t]\n' {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf
        bcftools sort -Oz -o {output.vcf} {output.vcf}.temp2.vcf
        tabix {output.vcf}
        rm {output.vcf}.temp*
    """

rule split_mini_vcf:
    input:
        vcf = rules.concat_quilt_vcf.output.vcf
    output:
        fv = f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/fv/quilt.chr{{chr}}.vcf.gz",
        mini = f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/mini/quilt.chr{{chr}}.vcf.gz",
        fv_unzip = temp(f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/fv/quilt.chr{{chr}}.vcf"),
        mini_unzip = temp(f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/mini/quilt.chr{{chr}}.vcf")
    resources:
        mem = '30G'
    params:
        fv = "data/sample_tsvs/chip_idt_names.tsv",
        fv_rename = "data/rename_tsvs/chip_idt_to_gam.tsv",
        mini = "data/sample_tsvs/mini_idt_names.tsv",
        mini_rename = "data/rename_tsvs/mini_idt_to_gam.tsv"
    shell: """ 
        mkdir -p results/mini_imputation/splited_vcfs/{PANEL_NAME}/fv/
        mkdir -p results/mini_imputation/splited_vcfs/{PANEL_NAME}/mini/

        bcftools view -S {params.fv} {input.vcf} | bcftools reheader -s {params.fv_rename} -o {output.fv_unzip}
        bcftools view -S {params.mini} {input.vcf} | bcftools reheader -s {params.mini_rename} -o {output.mini_unzip}

        bgzip {output.fv_unzip}
        bgzip {output.mini_unzip}

        touch {output.fv_unzip}
        touch {output.mini_unzip}
    """

rule filter_lc_info:
    input:
        lc_vcf = f"results/mini_imputation/splited_vcfs/{PANEL_NAME}/{{prep}}/quilt.chr{{chr}}.vcf.gz"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info/lc.chr{{chr}}.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        info = config['info_filter']
    shell: """
        mkdir -p results/wip_vcfs/{PANEL_NAME}/{wildcards.prep}/high_info/

        bcftools filter -i 'INFO_SCORE>{params.info}' -Oz -o {output.filtered_vcf} {input.lc_vcf}
    """

rule filter_lc_maf:
    input:
        lc_vcf = rules.filter_lc_info.output.filtered_vcf,
        af = "data/oneKG_MAFs/oneKG_AF_AFR_chr{chr}.txt"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info_high_af/lc.chr{{chr}}.vcf.gz"
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
             outdir="results/wip_vcfs/" + params.panel + "/" + wildcards.prep + "/high_info_high_af/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )

rule filter_lc_sites:
    input:
        vcf = rules.filter_lc_maf.output.filtered_vcf,
        sites = "data/chip/omni5m/omni5m_sites.tsv"
    output:
        filtered_vcf = f"results/wip_vcfs/{PANEL_NAME}/{{prep}}/high_info_high_af_chip_sites/lc.chr{{chr}}.vcf.gz"
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
             outdir="results/wip_vcfs/" + params.panel + "/" + wildcards.prep + "/high_info_high_af_chip_sites/",
             save_name="lc.chr" + str(wildcards.chr) + ".vcf.gz"
             )

pair = ['lc', 'hc']
axis = ['h', 'v']

imputation_dir = config['mini_imputation_dir']
lc_vcf_dir = config['mini_lc_vcf_dir']
hc_vcf_dir = config['mini_hc_vcf_dir']

def get_lc_vcf_dir(wildcards):
    ix = imputation_dir.index(wildcards.imp_dir)
    return lc_vcf_dir[ix]

def get_hc_vcf_dir(wildcards):
    ix = imputation_dir.index(wildcards.imp_dir)
    return hc_vcf_dir[ix]

rule copy_vcf_in_working_dir:
    input:
        lc_vcf_dir = get_lc_vcf_dir,
        hc_vcf_dir = get_hc_vcf_dir
    output:
        vcfs = temp(expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, allow_missing = True))
    resources:
        mem = '30G'
    threads: 4
    shell: """
        mkdir -p {wildcards.imp_dir}vcf/
        mkdir -p {wildcards.imp_dir}impacc/
        mkdir -p {wildcards.imp_dir}graphs/
        mkdir -p {wildcards.imp_dir}vcf/all_samples/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/all_samples/hc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_cc/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_cc/hc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_eth/lc_vcf/
        mkdir -p {wildcards.imp_dir}vcf/by_eth/hc_vcf/

        for i in {{1..22}}
        do
            cp {input.lc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/lc_vcf/lc.chr$i.vcf.gz
            cp {input.hc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/hc_vcf/hc.chr$i.vcf.gz
        done
    """

rule subset_lc_samples:
    input:
        quilt_vcf = '{imp_dir}vcf/all_samples/lc_vcf/lc.chr{chr}.vcf.gz',
        chip_vcf = '{imp_dir}vcf/all_samples/hc_vcf/hc.chr{chr}.vcf.gz'
    output:
        ss_lc_vcf = temp('{imp_dir}vcf/all_samples/lc_vcf/lc.subset.chr{chr}.vcf.gz'),
        ss_hc_vcf = temp('{imp_dir}vcf/all_samples/hc_vcf/hc.subset.chr{chr}.vcf.gz'),
        tmp_names = temp('{imp_dir}vcf/all_samples/samples.chr{chr}.tsv')
    resources:
        mem = '10G'
    run: 
        if (wildcards.imp_dir.split('/')[-3] == 'mini') and ('chip' in wildcards.imp_dir.split('/')[-2]):
            lc_names = lcwgsus.bcftools_get_samples(input.quilt_vcf)
            lcwgsus.save_lst(output.tmp_names, lc_names)
            shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_hc_vcf} {input.chip_vcf}")
            shell("cp {input.quilt_vcf} {output.ss_lc_vcf}")
        else:
            hc_names = lcwgsus.bcftools_get_samples(input.chip_vcf)
            lcwgsus.save_lst(output.tmp_names, hc_names)
            shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_lc_vcf} {input.quilt_vcf}")
            shell("cp {input.chip_vcf} {output.ss_hc_vcf}")

rule calculate_imputation_accuracy_all:
    input:
        quilt_vcf = rules.subset_lc_samples.output.ss_lc_vcf,
        chip_vcf = rules.subset_lc_samples.output.ss_hc_vcf,
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.tsv",
        h_impacc = "{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv",
        v_report = "{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.tsv",
        v_impacc = "{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv",
        lc_vcf = "{imp_dir}vcf/all_samples/filtered_vcfs/lc.chr{chr}.vcf.gz",
        hc_vcf = "{imp_dir}vcf/all_samples/filtered_vcfs/hc.chr{chr}.vcf.gz",
        af = "{imp_dir}vcf/all_samples/af/af.chr{chr}.tsv"
    resources:
        mem = '100G'
    threads: 16
    params:
        linker = config['sample_linker'],
        common_outdir = "{imp_dir}impacc/all_samples/",
        common_savedir = "{imp_dir}vcf/all_samples/",
        chrom = "{chr}"
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af_txt = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt, save_vcfs = True, lc_vcf_outdir = params.common_savedir + "filtered_vcfs/", hc_vcf_outdir = params.common_savedir + "filtered_vcfs/", lc_vcf_name = "lc.chr" + wildcards.chr + ".vcf.gz", hc_vcf_name = "hc.chr" + wildcards.chr + ".vcf.gz", af_outdir = params.common_savedir + "af/", af_name = "af.chr" + wildcards.chr + ".tsv")

        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = params.common_outdir + "by_variant/", 
                                                   save_name = 'chr' + wildcards.chr +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_variant/", 
                                           save_name = 'chr' + wildcards.chr +'.h.impacc.tsv')

        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = params.common_outdir + "by_sample/", 
                                           save_name = 'chr' + wildcards.chr +'.v.impacc.tsv')

        lcwgsus.rezip_vcf(output.lc_vcf)
        lcwgsus.rezip_vcf(output.hc_vcf)

rule aggregate_impaccs:
    input:
        h_impaccs = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        h_impacc = "{imp_dir}impacc/all_samples/by_variant/all.h.impacc.tsv",
        v_impacc = "{imp_dir}impacc/all_samples/by_sample/all.v.impacc.tsv"
    resources:
        mem = '30G'
    params:
        impacc_outdir = "{imp_dir}impacc/all_samples/"
    run:
        h_lst = input.h_impaccs
        outdir_h = params.impacc_outdir + "by_variant/"
        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs, save_file = True,
        outdir = outdir_h, save_name = 'all.h.impacc.tsv')

        v_lst = input.v_impaccs
        outdir_v = params.impacc_outdir + "by_sample/"
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs, save_file = True,
        outdir = outdir_v, save_name = 'all.v.impacc.tsv')

rule calculate_imputation_sumstat:
    input:
        h_impaccs = expand("{imp_dir}impacc/all_samples/by_variant/chr{chr}.h.impacc.tsv", chr = chromosome, allow_missing = True),
        v_impaccs = expand("{imp_dir}impacc/all_samples/by_sample/chr{chr}.v.impacc.tsv", chr = chromosome, allow_missing = True)
    output:
        sumstats = "{imp_dir}summary_metrics.tsv"
    resources:
        mem = '30G'
    params:
        axis = "v"
    run:
        lcwgsus.calculate_imputation_sumstats(wildcards.imp_dir, subset = False, axis = params.axis, save_file = True)