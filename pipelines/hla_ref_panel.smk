configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
sys.path.append("/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/")
import lcwgsus
from lcwgsus.variables import *
from hla_phase import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
samples_fv = read_tsv_as_lst("data/sample_tsvs/fv_idt_names.tsv")
chromosome = [i for i in range(1,23)]
QUILT_HOME = config["QUILT_HOME"]
NGEN = config["NGEN"]
RECOMB_POP = config["RECOMB_POP"]
studies = ['1KG', 'GAMCC']

'''
rule subset_vcf_to_chr6:
    input:
        vcf = "results/wip_vcfs/oneKG/vanilla/high_info_high_af_high_conf/lc.chr6.vcf.gz"
    output:
        tmp_vcf = temp("results/hla_tests/gamcc_vcf/fv.chr6.vcf.gz"),
        hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples"
    threads: 4
    resources: mem = '30G'
    params: 
        outdir = "results/hla_tests/gamcc_vcf/",
        fv = "data/sample_tsvs/fv_gm_names.tsv"
    shell: """
        mkdir -p {params.outdir}

        bcftools view -S {params.fv} {input.vcf}| \
        bcftools norm -m+ | \
        bcftools view -m2 -M2 -v snps | \
        
        bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert -h \
        {params.outdir}fv.chr6 {output.tmp_vcf}

        sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' {output.samples}
    """


rule prepare_hla_bamlist:
    input:
        bams = expand("data/bams/{id}.bam", id = samples_fv)
    output:
        bam_all = "results/hla_tests/bamlist.txt"
    localrule: True
    shell: """
        mkdir -p results/hla_tests/
        ls {input.bams} > {output.bam_all}
    """

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_versions = ['3390', '3570']

rule prepare_ref:
    input:
        hap = "results/hla_tests/gamcc_vcf/fv.chr6.hap.gz",
        legend = "results/hla_tests/gamcc_vcf/fv.chr6.legend.gz",
        samples = "results/hla_tests/gamcc_vcf/fv.chr6.samples",
        genetic_map = "data/imputation_accessories/maps/YRI-chr6-final.b38.txt"
    output:
        RData = "results/hla_tests/quilt.hrc.hla.all.haplotypes.RData"
    resources:
        mem = '30G'
    threads: 4
    params:
        outputdir = "results/hla_tests/prepared_ref/"
    shell: """
        mkdir -p {params.outputdir}
        R -e 'library("data.table"); library("QUILT");
        QUILT_prepare_reference(
        outputdir = "{params.outputdir}",
        nGen = {NGEN},
        chr = "chr6",
        regionStart = 25587319,
        regionEnd = 33629686,
        buffer = 500000,
        reference_haplotype_file = "{input.hap}",
        reference_legend_file = "{input.legend}",
        reference_sample_file = "{input.samples}",
        genetic_map_file = "{input.genetic_map}",
        reference_exclude_samplelist_file = "",
        output_file = "{output.RData}")'
    """

rule hla_imputation:
    input:
        bamlist = "results/hla/imputation/bamlists/bamlist{num}.txt",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_v{IPD_IMGT_version}/"
    output:
        imputed = "results/hla/imputation/QUILT_HLA_result_v{IPD_IMGT_version}/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '50G'
    threads: 6
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_v{IPD_IMGT_version}/genes{wildcards.num}/{wildcards.hla_gene}/

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_v{wildcards.IPD_IMGT_version}/genes{wildcards.num}/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={input.ref_dir} \
        --quilt_hla_haplotype_panelfile={input.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """
'''

rule phase_GAMCC_alleles:
    input:
        phased_vcf_file = "results/imputation/vcfs/malariaGen_v1_b38/quilt.chr6.vcf.gz"
    output:
        html = "results/hla_tests/phasing/html/GAMCC-{gene}.html",
        phase_df = "results/hla_tests/phasing/phased_dfs/GAMCC-{gene}.tsv"
    resources:
        mem = '30G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv',
        reference_allele_file = '/well/band/users/rbx225/recyclable_files/hla/b38_reference_alleles.tsv',
        gm_tsv_names = 'data/sample_tsvs/fv_gm_names.tsv'
    run:
        gene = wildcards.gene
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')

        gamcc_hla = lcwgsus.read_hla_direct_sequencing(retain = 'fv', unique_two_field = False)
        gamcc_hla = gamcc_hla[['SampleID', 'Locus', 'Two field1', 'Two field2']].reset_index(drop = True)

        colnames = ['Sample ID'] + [label for g in HLA_GENES for label in [f'HLA-{g} 1', f'HLA-{g} 2']]
        hlatypes = pd.DataFrame(columns = colnames)
        for s in gamcc_hla['SampleID'].unique():
            tmp = gamcc_hla[gamcc_hla['SampleID'] == s]
            row = [s] + tmp[['Two field1', 'Two field2']].values.ravel().tolist()
            hlatypes.loc[len(hlatypes)] = row
    
        reference_allele_ary = np.array(lcwgsus.read_tsv_as_lst(params.reference_allele_file))

        subset_vcf_samples = lcwgsus.read_tsv_as_lst(params.gm_tsv_names)
        subset_vcf_samples = ','.join(subset_vcf_samples)

        sample_linker = pd.read_csv(SAMPLE_LINKER_FILE)
        sample_linker = {k:v for k, v in zip(sample_linker['Sample_Name'], sample_linker['Chip_Name'])}

        return_dict = phase_hla_on_haplotypes(gene = gene, 
                            ipd_gen_file_dir = params.ipd_gen_file_dir, 
                            hla_gene_information = hla_gene_information, 
                            hlatypes = hlatypes,
                            phased_vcf = input.phased_vcf_file, 
                            reference_allele_ary = reference_allele_ary, 
                            strict_snp_filter = False,
                            read_from_QUILT = True, 
                            subset_vcf_samples = subset_vcf_samples,
                            sample_linker = sample_linker)
        hlatypes = return_dict['hlatypes']

        individual = 'GAM951638'
        ix = hlatypes.index[hlatypes['Sample ID'] == individual][0]
        display_indices = np.arange(10)

        res = visualise_phase(gene, ix, hlatypes, return_dict, both_het = True)
        compare_phase(display_indices, res, save_html = True, save_name = output.html)
        df = return_dict['phase_df'][['Sample', 'allele1', 'allele2']]
        df.columns = ['Sample ID', f'HLA-{gene} 1', f'HLA-{gene} 2']
        df.to_csv(output.phase_df, sep = '\t', index = False, header = True)

rule phase_1KG_alleles:
    input:
        samples_file = "results/hla/imputation/ref_panel/auxiliary_files/oneKG.samples",
        hlatypes_file = "results/hla/imputation/ref_panel/auxiliary_files/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        phased_vcf_file = "/well/band/users/rbx225/recyclable_files/ref_panels/oneKG/oneKG.chr6.vcf.gz"
    output:
        html = "results/hla_tests/phasing/html/1KG-{gene}.html",
        phase_df = "results/hla_tests/phasing/phased_dfs/1KG-{gene}.tsv"
    resources:
        mem = '30G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv',
        reference_allele_file = '/well/band/users/rbx225/recyclable_files/hla/b38_reference_alleles.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        ref_samples = pd.read_csv(input.samples_file, sep = ' ')
        hlatypes = pd.read_csv(input.hlatypes_file, sep = '\t')

        ref_samples_removed = ref_samples[~ref_samples['SAMPLE'].isin(hlatypes['Sample ID'].tolist())]
        samples_to_remove = ref_samples_removed['SAMPLE'].tolist()
        hlatypes = hlatypes[~hlatypes['Sample ID'].isin(samples_to_remove)].sort_values(by = 'Sample ID').reset_index(drop = True)

        reference_allele_ary = np.array(lcwgsus.read_tsv_as_lst(params.reference_allele_file))

        return_dict = phase_hla_on_haplotypes(gene = wildcards.gene, 
                                    ipd_gen_file_dir = params.ipd_gen_file_dir, 
                                    hla_gene_information = hla_gene_information,
                                    hlatypes = hlatypes,
                                    phased_vcf = input.phased_vcf_file, 
                                    reference_allele_ary = reference_allele_ary, 
                                    strict_snp_filter = True,
                                    read_from_QUILT = False, 
                                    subset_vcf_samples = None,
                                    sample_linker = None)
        hlatypes = return_dict['hlatypes']

        individual = 'NA12878'
        ix = hlatypes.index[hlatypes['Sample ID'] == individual][0]
        display_indices = np.arange(10)

        res = visualise_phase(gene, ix, hlatypes, return_dict, both_het = True)
        compare_phase(display_indices, res, save_html = True, save_name = output.html)
        df = return_dict['phase_df'][['Sample', 'allele1', 'allele2']]
        df.columns = ['Sample ID', f'HLA-{gene} 1', f'HLA-{gene} 2']
        df.to_csv(output.phase_df, sep = '\t', index = False, header = True)


to_merge = ['gamcc', 'oneKG']

rule pre_prepare_merge_GAMCC_vcf:
    input:
        vcf = "results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.dose.vcf.gz"
    output:
        fv_vcf = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf.gz",
        tmp_vcf = temp("results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr{chr}.vcf")
    resources: mem = '30G'
    threads: 4
    params: 
        outdir = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/",
        fv_gm_names = "data/sample_tsvs/fv_gm_names.tsv",
        gm_to_gam = "data/rename_tsvs/fv_gm_to_gam.ssv"
    shell: """
        mkdir -p {params.outdir}

        bcftools view -S {params.fv_gm_names} {input.vcf} | \
        bcftools reheader -s {params.gm_to_gam} -o {output.tmp_vcf}

        bgzip {output.tmp_vcf}
    """

def get_vcf_from_to_merge(wildcards):
    if wildcards.to_merge == 'oneKG':
        return "data/ref_panel/oneKG/oneKG.chr" + wildcards.chr + ".vcf.gz"
    elif wildcards.to_merge == 'gamcc':
        return "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr" + wildcards.chr + ".vcf.gz"
    else:
        return ""

rule prepare_merge_1KG_GAMCC_vcf:
    input:
        vcf = get_vcf_from_to_merge
    output:
        haps = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.hap"),
        legend = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.legend")
    resources: mem = '30G'
    threads: 4
    params: outdir = "results/hla_ref_panel/oneKG_mGenv1/tmp/"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | \
        bcftools view -m2 -M2 -v snps | \
        bcftools sort | \
        bcftools convert -h {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}

        gunzip {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.hap.gz
        gunzip {params.outdir}{wildcards.to_merge}.chr{wildcards.chr}.legend.gz
    """

rule prepare_merge_1KG_GAMCC_sample:
    input:
        mg = "results/hla_ref_panel/oneKG_mGenv1/fv_gamcc_vcf/gamcc.chr22.vcf.gz",
        oneKG = "data/ref_panel/oneKG/oneKG.chr22.vcf.gz"
    output:
        tmp_sample = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/tmp.sample"),
        sample = temp("results/hla_ref_panel/oneKG_mGenv1/tmp/merged.samples")
    resources: mem = '30G'
    threads: 1
    params: outdir = "results/hla_ref_panel/oneKG_mGenv1/tmp/"
    shell: """
        mkdir -p {params.outdir}

        echo "sample population group sex" >> {output.sample}

        bcftools query -l {input.oneKG} >> {output.tmp_sample}
        bcftools query -l {input.mg} >> {output.tmp_sample}

        for i in $(cat {output.tmp_sample}); do
            echo $i $i $i 2 >> {output.sample}
        done
    """

# Maybe instead of per_chunk use only the HLA region?
rule merge_1KG_GAMCC_per_chunk:
    input:
        haps = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.hap", to_merge = to_merge, allow_missing = True),
        legends = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr{chr}.legend", to_merge = to_merge, allow_missing = True),
        gen_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr{{chr}}-final.b38.txt",
        sample = rules.prepare_merge_1KG_GAMCC_sample.output.sample
    output:
        haps = temp("results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.hap"),
        legend = temp("results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.legend"),
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}.vcf.gz",
    resources: mem = '70G'
    threads: 4
    params: 
        impute2 = tools['impute2'],
        outdir = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/",
        output_prefix = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr{chr}.{regionStart}.{regionEnd}",
        mGen_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr{chr}.hap",
        mGen_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr{chr}.legend",
        oneKG_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr{chr}.hap",
        oneKG_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr{chr}.legend",
    shell: """
        mkdir -p {params.outdir}

       {params.impute2} \
        -merge_ref_panels_output_ref {params.output_prefix}.tmp \
        -m {input.gen_map} \
        -h {params.oneKG_haps} \
           {params.mGen_haps} \
        -l {params.oneKG_legend} \
           {params.mGen_legend} \
        -int {wildcards.regionStart} {wildcards.regionEnd} \
        -k_hap 2018 420 \
        -Ne 20000

        awk -F ' ' 'NR==1 {{print; next}} {{$1 = "chr{wildcards.chr}:"$2"_"$3"_"$4; print $0}}' \
        {params.output_prefix}.tmp.legend > {output.legend}
        mv {params.output_prefix}.tmp.hap {output.haps}

        bgzip -f {output.legend}
        bgzip -f {output.haps}
        touch {output.legend}
        touch {output.haps}

        cp {input.sample} {params.output_prefix}.samples

        bcftools convert -H {params.output_prefix} | bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

def get_input_vcfs_as_list(wildcards):
    return(mGen_chunk_vcf_dict[str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, mGen_chunk_vcf_dict[str(wildcards.chr)])))

rule merge_1KG_GAMCC_chunks:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/oneKG_GAMCC.chr{chr}.vcf.gz",
    resources: mem = '30G'
    threads: 4
    params: 
        input_string = get_input_vcfs_as_string,
    shell: """
        bcftools concat --ligate-force -Oz -o {output.vcf} {params.input_string}
        tabix {output.vcf}
    """

hla_ref_panel_start = config["hla_b38_start_extra"]
hla_ref_panel_end = config["hla_b38_end_extra"]

rule merge_1KG_GAMCC_hla_only:
    input:
        haps = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr6.hap", to_merge = to_merge),
        legends = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr6.legend", to_merge = to_merge),
        gen_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr{{chr}}-final.b38.txt",
        sample = rules.prepare_merge_1KG_GAMCC_sample.output.sample
    output:
        haps = temp(f"results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr6.{hla_ref_panel_start}.{hla_ref_panel_end}.hap"),
        legend = temp(f"results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr6.{hla_ref_panel_start}.{hla_ref_panel_end}.legend"),
        vcf = "results/hla_ref_panel/oneKG_mGenv1/merged/hla/chr6.hla.vcf.gz",
    resources: mem = '60G'
    threads: 4
    params: 
        impute2 = tools['impute2'],
        outdir = "results/hla_ref_panel/oneKG_mGenv1/merged/regions/",
        output_prefix = f"results/hla_ref_panel/oneKG_mGenv1/merged/regions/chr6.{hla_ref_panel_start}.{hla_ref_panel_end}",
        mGen_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr6.hap",
        mGen_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/gamcc.chr6.legend",
        oneKG_haps = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr6.hap",
        oneKG_legend = "results/hla_ref_panel/oneKG_mGenv1/tmp/oneKG.chr6.legend",
    shell: """
        mkdir -p {params.outdir}

       {params.impute2} \
        -merge_ref_panels_output_ref {params.output_prefix}.tmp \
        -m {input.gen_map} \
        -h {params.oneKG_haps} \
           {params.mGen_haps} \
        -l {params.oneKG_legend} \
           {params.mGen_legend} \
        -int {hla_ref_panel_start} {hla_ref_panel_end} \
        -k_hap 2018 420 \
        -Ne 20000

        awk -F ' ' 'NR==1 {{print; next}} {{$1 = "chr6:"$2"_"$3"_"$4; print $0}}' \
        {params.output_prefix}.tmp.legend > {output.legend}
        mv {params.output_prefix}.tmp.hap {output.haps}

        bgzip -f {output.legend}
        bgzip -f {output.haps}
        touch {output.legend}
        touch {output.haps}

        cp {input.sample} {params.output_prefix}.samples

        bcftools convert -H {params.output_prefix} | bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """