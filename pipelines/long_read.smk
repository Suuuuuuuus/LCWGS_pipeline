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

read_lengths = ['1kb', '2kb', '5kb', '10kb', '20kb']
# read_lengths = ['1kb']
haplotypes = ['mat', 'pat']
# haplotypes = ['mat']
method = 'CCS'
coverage = '0.5'

def get_num_mean_length(wildcards):
    return int(wildcards.rl[:-2])*1000

rule simulate_reads:
    input:
        fasta = "data/lr_fasta/HG02886.{hap}.fa"
    output:
        tmp = temp("data/lr_simulations/{rl}/{hap}.{rl}.fastq"),
        fastq = "data/lr_simulations/{rl}/{hap}.{rl}.fastq.gz"
    resources:
        mem = '30G'
    threads: 4
    params:
        model = "data/lr_models/model_qc_" + method.lower(),
        mean_length = get_num_mean_length,
        max_length = 50000,
        outdir = "data/lr_simulations/{rl}/",
        prefix = "data/lr_simulations/{rl}/tmp.{hap}.{rl}"
    shell: """
        mkdir -p {params.outdir}

        pbsim --data-type {method} \
        --depth {coverage} \
        --model_qc {params.model} \
        --length-mean {params.mean_length} \
        --length-max {params.max_length} \
        --prefix {params.prefix} {input.fasta}

        cat {params.prefix}*.fastq > {output.tmp}

        rm {params.prefix}*.fastq
        rm {params.prefix}*.maf
        rm {params.prefix}*.ref

        gzip {output.tmp}
        touch {output.tmp}
    """

rule lr_alignment:
    input:
        fastq = rules.simulate_reads.output.fastq,
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if config['concatenate'] else config["ref38"]
    output:
        sam = temp("data/lr_bams/tmp/{hap}.{rl}.sam"),
        bam = temp("data/lr_bams/tmp/{hap}.{rl}.bam")
    resources:
        mem = '30G'
    threads: 6
    shell: """
        mkdir -p data/lr_bams/tmp/

        minimap2 -t {threads} -ax map-hifi {input.reference} {input.fastq} > {output.sam}
        samtools view -bS {output.sam} > {output.bam}
    """

#         pbmm2 align {input.reference} {input.fastq} {output.bam} \
        # --sort -j 4 -J 2

rule lr_clean_bam:
    input:
        bams = expand("data/lr_bams/tmp/{hap}.{rl}.bam", hap = haplotypes, allow_missing = True)
    output:
        bam = "data/lr_bams/{rl}.bam",
        bai = "data/lr_bams/{rl}.bam.bai",
        tmp1 = temp("data/lr_bams/{rl}.tmp1.bam"),
        metric = temp("data/lr_bams/{rl}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/lr_bams/tmp/{rl}/",
        sample = "{rl}"
    shell: """
        mkdir -p {params.tmpdir}

        samtools cat -o {output.tmp1} {input.bams}
        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        samtools index {output.bam}

        picard AddOrReplaceReadGroups \
        -VERBOSITY ERROR \
        -I {output.bam} \
        -O {output.tmp1} \
        -RGLB OGC \
        -RGPL PacBio \
        -RGPU unknown \
        -RGSM {params.sample}

        java -Xmx40G -Xms20G -jar /well/band/users/rbx225/conda/skylake/envs/sus/share/picard-slim-2.27.4-0/picard.jar \
        FixMateInformation -I {output.tmp1}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        java -Xmx40G -Xms20G -jar /well/band/users/rbx225/conda/skylake/envs/sus/share/picard-slim-2.27.4-0/picard.jar \
        MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

QUILT_HOME = config["QUILT_HOME"]
lr_analysis_dir = config["lr_analysis_dir"]
RECOMB_POP=config["RECOMB_POP"]
NGEN=config["NGEN"]
WINDOWSIZE=config["WINDOWSIZE"]
BUFFER=config["BUFFER"]
PANEL_NAME=config["PANEL_NAME"]

rule prepare_bamlist:
    input:
        bams = expand("data/lr_bams/{rl}.bam", rl = read_lengths)
    output:
        bamlist = "results/lr_imputation/bamlist.txt"
    shell: """
        mkdir -p {lr_analysis_dir}

        ls data/lr_bams/*.bam > {output.bamlist}
    """

rule convert_recomb:
    input:
        f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz"
    output:
        f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        threads = 1
    wildcard_constraints:
        chr='\d{1,2}'
    shell: """
        R -f {QUILT_HOME}scripts/make_b38_recomb_map.R \
        --args {lr_analysis_dir} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{PANEL_NAME}/{PANEL_NAME}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/lr_imputation/refs/tmp.{PANEL_NAME}.chr{{chr}}.vcf.gz"),
        hap = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        samples = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    params:
        panel = PANEL_NAME,
        threads=1
    shell: """
        mkdir -p results/lr_imputation/refs/
        bcftools norm -m+ {input.vcf} | bcftools view -Oz -o {output.tmp_vcf} -m2 -M2 -v snps

        tabix {output.tmp_vcf}

        bcftools convert --haplegendsample \
        results/lr_imputation/refs/{params.panel}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz", chr = chromosome),
        code = "scripts/determine_chunks.R"
    output:
        json = "results/lr_imputation/regions.json"
    resources: mem = '10G'
    shell: """
        Rscript {input.code} {lr_analysis_dir:q} {WINDOWSIZE} {BUFFER} {PANEL_NAME:q}
    """

rule prepare_ref:
    input:
        json = "results/lr_imputation/regions.json",
        hap = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.hap.gz",
        legend = f"results/lr_imputation/refs/{PANEL_NAME}.chr{{chr}}.legend.gz",
        recomb = f"results/lr_imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/lr_imputation/refs/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem_mb = 30000
    params:
        threads = 8
    shell: """
        mkdir -p results/lr_imputation/refs/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/lr_imputation/refs/RData/other/", \
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
        bamlist = "results/lr_imputation/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = f"results/lr_imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
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
        mkdir -p "results/lr_imputation/vcfs/{params.panel}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/lr_imputation/refs/RData/other/", \
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

file="results/lr_imputation/regions.json"
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
        file="results/lr_imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)]=per_chr_vcfs
    final_vcfs.append("results/lr_imputation/vcfs/" + PANEL_NAME + "/quilt.chr" + str(chr) + ".vcf.gz")

def get_input_vcfs_as_list(wildcards):
    return(vcfs_to_concat[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, vcfs_to_concat[str(wildcards.chr)])))

rule concat_quilt_vcf:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/lr_imputation/vcfs/{PANEL_NAME}/quilt.chr{{chr}}.vcf.gz"
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

rule prepare_lr_vcf:
    input:
        vcf = f"data/ref_panel/{PANEL_NAME}/oneKG.chr{{chr}}.vcf.gz"
    output:
        one = temp("results/lr_imputation/truth/1kb.chr{chr}.vcf.gz"),
        two = temp("results/lr_imputation/truth/2kb.chr{chr}.vcf.gz"),
        five = temp("results/lr_imputation/truth/5kb.chr{chr}.vcf.gz"),
        ten = temp("results/lr_imputation/truth/10kb.chr{chr}.vcf.gz"),
        twenty = temp("results/lr_imputation/truth/20kb.chr{chr}.vcf.gz"),
        truth = "results/lr_imputation/truth/long_read_truth.chr{chr}.vcf.gz",
        rename = temp("results/lr_imputation/truth/name.chr{chr}.txt")
    resources: mem = '10G'
    params:
        sample = 'HG02886'
    shell: """
        mkdir -p results/lr_imputation/truth/
        
        length=("10kb" "1kb" "20kb" "2kb" "5kb")

        for l in "${{length[@]}}"
        do
            echo $l > {output.rename}
            bcftools view -s {params.sample} {input.vcf} | \
            bcftools reheader -s {output.rename} -o results/lr_imputation/truth/$l.chr{wildcards.chr}.vcf
            bgzip results/lr_imputation/truth/$l.chr{wildcards.chr}.vcf
            tabix results/lr_imputation/truth/$l.chr{wildcards.chr}.vcf.gz
        done

        bcftools merge -Oz -o {output.truth} {output.one} {output.two} {output.five} {output.ten} {output.twenty}
    """

rule copy_vcf_in_working_dir_hack:
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
            cp /well/band/users/rbx225/GAMCC/results/lr_imputation/vcfs/oneKG/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/lc_vcf/lc.chr$i.vcf.gz
            cp /well/band/users/rbx225/GAMCC/results/lr_imputation/truth/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/hc_vcf/hc.chr$i.vcf.gz
        done
    """

pair = ['lc', 'hc']
axis = ['h', 'v']

# imputation_dir = list(config['imputation_dir'][-1])
# lc_vcf_dir = list(config['lc_vcf_dir'][-1])
# hc_vcf_dir = list(config['hc_vcf_dir'][-1])

# def get_lc_vcf_dir(wildcards):
#     ix = imputation_dir.index(wildcards.imp_dir)
#     return lc_vcf_dir[ix]

# def get_hc_vcf_dir(wildcards):
#     ix = imputation_dir.index(wildcards.imp_dir)
#     return hc_vcf_dir[ix]

# rule copy_vcf_in_working_dir:
#     input:
#         lc_vcf_dir = get_lc_vcf_dir,
#         hc_vcf_dir = get_hc_vcf_dir
#     output:
#         vcfs = temp(expand('{imp_dir}vcf/all_samples/{pair}_vcf/{pair}.chr{chr}.vcf.gz', chr = chromosome, pair = pair, allow_missing = True))
#     resources:
#         mem = '30G'
#     threads: 4
#     shell: """
#         mkdir -p {wildcards.imp_dir}vcf/
#         mkdir -p {wildcards.imp_dir}impacc/
#         mkdir -p {wildcards.imp_dir}graphs/
#         mkdir -p {wildcards.imp_dir}vcf/all_samples/lc_vcf/
#         mkdir -p {wildcards.imp_dir}vcf/all_samples/hc_vcf/
#         mkdir -p {wildcards.imp_dir}vcf/by_cc/lc_vcf/
#         mkdir -p {wildcards.imp_dir}vcf/by_cc/hc_vcf/
#         mkdir -p {wildcards.imp_dir}vcf/by_eth/lc_vcf/
#         mkdir -p {wildcards.imp_dir}vcf/by_eth/hc_vcf/

#         for i in {{1..22}}
#         do
#             cp {input.lc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/lc_vcf/lc.chr$i.vcf.gz
#             cp {input.hc_vcf_dir}/*chr$i.*.gz {wildcards.imp_dir}vcf/all_samples/hc_vcf/hc.chr$i.vcf.gz
#         done
#     """

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
        if (wildcards.imp_dir.split('/')[-3] == 'mini') and ('lc_chip' in wildcards.imp_dir.split('/')[-2]):
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