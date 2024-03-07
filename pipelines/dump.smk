configfile: "pipelines/config.json"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/band/users/rbx225/software/lcwgsus/")
import lcwgsus

config['samples'] = pd.read_table("samples.tsv", header = None, names = ['Code'])
ids_1x_all = list(config['samples']['Code'].values)
chromosome = [i for i in range(1,23)]

subsample_coverage = config['subsample_depth']
rm_bed_regions = config['rm_bed_regions']
bed_regions = config['bed_regions']

rule multiqc_hc:
    input:
        html1 = expand("results/fastqc/{id}_1_fastqc.html",id = samples_hc),
        html2 = expand("results/fastqc/{id}_2_fastqc.html",id = samples_hc),
        zip = expand("results/fastqc/{id}_{read}_fastqc.zip",id = samples_hc, read = ['1', '2'])
    output:
        html = "results/fastqc/multiqc_hc/multiqc_report.html",
        directory("results/fastqc/multiqc_hc")
    threads: 1
    resources:
        mem = '10G'
    params:
        outdir = "results/fastqc/multiqc_hc"
    shell:
        "multiqc {input.zip} --interactive -o {params.outdir}"

# The following rule was used to chunk-up bam files into chromosomes, but currently deprecated as GATK is very strict about bam files and cannot accept singleton reads
sample_linker = pd.read_table(config['sample_linker'], sep = ',')
ids_1x_all = list(sample_linker['Seq_Name'].values) # to be deprecated
test_hc = ids_1x_all[:2]
test_hc_dict = read_tsv_as_dict(test_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")
test_hc_all_chunks = []
for value_list in test_hc_dict.values():
    test_hc_all_chunks.extend(value_list)
rule split_bams:
    input:
        bam = "data/bams/{id}.bam",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta"
    output:
        bam_chunk = temp("data/chunk_bams/tmp/tmp/{id}/{id}.chr{chr}.bam"),
        tmp1 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.chr{chr}.tmp1.bam"),
        tmp2 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.chr{chr}.tmp2.bam")
#        reference = rules.index_reference.input.reference
    threads: 4
    resources: mem = '10G'
    params:
        chr_str = "chr{chr}",
        verbosity = "ERROR",
        sample = "{id}".split("_")[0]
    shell: """
        mkdir -p data/chunk_bams/tmp/tmp/{wildcards.id}/
        samtools view -h {input.bam} {params.chr_str} | \
        samtools sort -o {output.tmp1} -

        samtools index {output.tmp1}

        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {output.tmp1} \
        -O {output.tmp2} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {params.sample}

        picard FixMateInformation -I {output.tmp2}

        samtools collate -Oun128 {output.tmp2} | samtools fastq -OT RG,BC - | \
        bwa mem -pt4 -CH <(samtools view -H {output.tmp2} | grep ^@RG) {input.reference} - | \
        samtools sort -@4 -m4g -o {output.bam_chunk} -
    """

rule calculate_per_bin_coverage_1x:
    input:
        script = "scripts/calculate_per_bin_coverage.py",
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph.bed"
    output:
        coordinates = "results/coverage/per_bin_coverage/1x/{id}_chr{chr}_coordinate.txt",
        bases = "results/coverage/per_bin_coverage/1x/{id}_chr{chr}_base.txt"
    params:
        indir = "results/coverage/bedgraphs/",
        outdir = "results/coverage/per_bin_coverage/1x/",
        bin_size = 10000
    threads: 8
    resources: mem_mb = 50000
    shell: """
        python {input.script} {wildcards.id} {wildcards.chr} {params.indir} {params.outdir} {params.bin_size}
    """

rule calculate_per_bin_coverage_20x:
    input:
        script = "scripts/calculate_per_bin_coverage.py"
    output:
        coordinates = "results/coverage/per_bin_coverage/20x/{id_20x}_chr{chr}_coordinate.txt",
        bases = "results/coverage/per_bin_coverage/20x/{id_20x}_chr{chr}_base.txt"
    params:
        indir = "results/coverage/bedgraphs/",
        outdir = "results/coverage/per_bin_coverage/20x/",
        bin_size = 10000
    # threads: 16
    # resources: mem_mb = 40000
    shell: """
        python {input.script} {wildcards.id_20x} {wildcards.chr} {params.indir} {params.outdir} {params.bin_size}
    """

rule samtools_coverage:
    input:
        bam = "data/bams/{id}.bam"
    output:
        per_chromosome_coverage = "results/coverage/per_chromosome_coverage/{id}_per_chromosome_coverage.txt"
    shell: """
        samtools coverage {input.bam} | sed -n '2,23p' > {output.per_chromosome_coverage}
    """

rule calculate_uncoverage_rate:
    input:
        per_chromosome_coverage = rules.samtools_coverage.output.per_chromosome_coverage
    output:
        uncoverage_rate = temp("results/coverage/per_chromosome_coverage/{id}_uncoverage_rate.txt")
    shell: """
	total=$(cut -f3 {input.per_chromosome_coverage} | paste -sd+ | bc)
        covered=$(cut -f5 {input.per_chromosome_coverage} | paste -sd+ | bc)
        result=$(echo "scale=4; (1-$covered/$total)" | bc)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
    """

rule extract_fastqc_dup_rate:
    input:
        zip1 = "results/fastqc/{id}_1_fastqc.zip",
        zip2 = "results/fastqc/{id}_2_fastqc.zip"
    output:
        txt = temp("results/fastqc/tmp/{id}.txt")
    params:
        tmpdir = "results/fastqc/"
    shell: """
        unzip "{params.tmpdir}{wildcards.id}_1_fastqc.zip" -d {params.tmpdir}
        unzip "{params.tmpdir}{wildcards.id}_2_fastqc.zip" -d {params.tmpdir}
        dedup1=$(grep "Total Deduplicated Percentage" "{params.tmpdir}{wildcards.id}_1_fastqc/fastqc_data.txt" | cut -f 2)
        dedup2=$(grep "Total Deduplicated Percentage" "{params.tmpdir}{wildcards.id}_2_fastqc/fastqc_data.txt" | cut -f 2)
        echo "{wildcards.id} $dedup1 $dedup2" > {output.txt}
        rm -r {params.tmpdir}*{wildcards.id}*fastqc
    """

rule aggregate_fastqc_dup_rate:
    input:
        files = expand("results/fastqc/tmp/{id}.txt", id = ids_1x_all)
    output:
        fastqc_dup_rate = temp("results/fastqc/duplication_rate_fastqc_tmp.txt")
    params:
        tmpdir = "results/fastqc/tmp/"
    shell: """
        cat {input.files} >> {output.fastqc_dup_rate}
        rm -r {params.tmpdir}
    """

rule calculate_fastqc_dup_rate:
    input:
        fastqc_dup_rate = rules.aggregate_fastqc_dup_rate.output.fastqc_dup_rate
    output:
        fastqc_dup_rate_res = "results/fastqc/duplication_rate_fastqc.txt"
    shell: """
        while IFS=$' ' read -r col1 col2 col3; do
            result=$(echo "scale=4; (1-($col3+$col2)/200)" | bc)
            formatted_result=$(printf "%06.4f" $result)
            echo -e "$col1\t$formatted_result" >> {output.fastqc_dup_rate_res}
        done < {input.fastqc_dup_rate}
    """

rule quilt_info:
    input:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    output:
        vcf = f"results/imputation/vcfs/{PANEL_NAME}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz.output.RData"
    params:
        threads = 1
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        R -f ${{QUILT_WRAP_HOME}}info.R --args {output.vcf}
    """

rule aggregate_results:
    input:
        script = "scripts/lcwgs_result_wrap_up.py",
        fastqc_dup_rate = "results/fastqc/duplication_rate_fastqc.txt",
        uncoverage_rate = "results/coverage/per_chromosome_coverage/uncoverage_rate.txt",
        samtools_dup_rate = "results/dup_rate/duplication_rate_samtools.txt",
        kmer_accuracy1 = "results/kmer/kmer_accuracy_read1.txt",
        kmer_accuracy2 = "results/kmer/kmer_accuracy_read2.txt",
        coverage = "results/coverage/per_sample_coverage.txt",
        fragment_size = "results/fragment_size/fragment_size.txt",
        proportion_ss_fragment_size = "results/fragment_size/proportion_ss_fragment_size.txt",
        proportion_fragment_size = "results/fragment_size/proportion_fragment_size.txt",
        fragment_overlap = "results/fragment_size/fragment_overlap.txt"
    output:
        result = "results/lcwgs_results.csv"
    shell: """
        python {input.script}
    """

# Generate bed chunk files from UCSC genome sizes
rule get_bam_bed_chunks: # Need to resolve the dependence with split_bams
    input:
        bed = "data/bedgraph/GRCh38.autosomes.bed"
    output:
        bed_chunks_samtools = "data/bedgraph/bam_chunks.bed", # samtools use chr1:1-2
        bed_chunks_names = "data/bedgraph/bam_chunks_names.bed" # names use chr1-1-2 to avoid naming errors
    threads: 1
    params:
        bam_chunk_size = config["bam_chunk_size"]
    shell: """
        bedtools makewindows -b {input.bed} -w {params.bam_chunk_size} | awk '{{print $1":"$2"-"$3}}' > {output.bed_chunks_samtools}
        sed 's/:/-/g' {output.bed_chunks_samtools} > {output.bed_chunks_names}
    """

rule split_bams:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bam_chunk = temp("data/chunk_bams/tmp/tmp/{id}/{id}.chr{chr}.bam")
    threads: 1
    resources: mem = '10G'
    params:
        chr_str = "chr{chr}"
    shell: """
        mkdir -p data/chunk_bams/tmp/tmp/{wildcards.id}/
        samtools view -h {input.bam} {params.chr_str} | \
        samtools sort -n - | \
        samtools fixmate -m - - -u | \
        samtools sort - -u | \
        samtools markdup - {output.bam_chunk}
    """

rule merge_bam:
    input:
        bams = lambda wildcards: nest[str(wildcards.hc)][str(wildcards.chr)]
    output:
        bam = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"), ### Now the wildcards are messed up to avoid intermediate files... Need to come back later
        bai = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai")
    threads: 2
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/chunk_bams/tmp/{wildcards.hc}/
        
        samtools merge -u - {input.bams} | \
        samtools sort -n - | \
        samtools fixmate -m - - -u | \
        samtools sort - -u | \
        samtools markdup - {output.bam}

        samtools index {output.bam}
    """

test_hc_dict = read_tsv_as_dict(test_hc, "data/file_lsts/hc_fastq_split/", "_split.tsv")
nest = {}
for i in test_hc:
    nest[i] = {}
    for j in chromosome:
        nest[str(i)][str(j)] = ["data/chunk_bams/tmp/tmp/" + k + "/" + k + ".chr" + str(j) + ".bam" for k in test_hc_dict[str(i)]]
rule merge_splited_bam:
    input:
        bams = lambda wildcards: nest[str(wildcards.hc)][str(wildcards.chr)]
    output:
        # bam = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"), ### Now the wildcards are messed up to avoid intermediate files... Need to come back later
        # bai = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai")
        bam = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"
#        bai = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam.bai"
    threads: 2
    resources:
        mem = '50G'
    shell: """
        mkdir -p data/chunk_bams/tmp/{wildcards.hc}/
        samtools cat -o {output.bam} {input.bams}
    """

rule rmdup_split:
    input:
        bam_chunk = "data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.bam"
    output:
        tmp1 = temp("data/chunk_bams/{hc}/{hc}.chr{chr}.tmp1.bam"),
        tmp2 = temp("data/chunk_bams/{hc}/{hc}.chr{chr}.tmp2.bam"),
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        dedup_bai_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam.bai",
        metric = temp("data/chunk_bams/tmp/{hc}/{hc}.chr{chr}.metrics.txt")
    threads: 2
    resources: mem = '50G'
    shell: """
        samtools sort -o {output.tmp1} {input.bam_chunk}

        picard MarkDuplicates \
        -I {output.tmp1} \
        -O {output.tmp2} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -o {output.dedup_bam_chunk} {output.tmp2}

        samtools index {output.dedup_bam_chunk}
    """

### variant calling ###

rule GATK_add_readgroup:
    input:
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
    output:
    
    params:
        verbosity = "ERROR"
    shell: """
        picard AddOrReplaceReadGroups \
        -VERBOSITY {params.verbosity} \
        -I {output.bam_marked_dups} \
        -O {output.bam_read_group} \
        -RGLB OGC \
        -RGPL ILLUMINA \
        -RGPU unknown \
        -RGSM {wildcards.hc}
    """

rule get_bqsr_report:
    input:
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        bqsr_report = "results/call/BQSR/BQSR_reports/{hc}.chr{chr}.BQSR.report"
    params:
        bqsr_known_sites = config["bqsr_known_sites"]
    run:
        cmd = ""
        for file in params.bqsr_known_sites:
            cmd = cmd + "--known-sites " + file + " "
        shell("""
            gatk --java-options "-Xmx8G" BaseRecalibrator \
            -I {dedup_bam_chunk} \
            -R {ref} \
            -O {report} \
            {cmd}
        """.format(cmd = cmd, dedup_bam_chunk = input.dedup_bam_chunk, ref = input.reference, report = output.bqsr_report))

rule apply_bqsr:
    input:
        dedup_bam_chunk = "data/chunk_bams/{hc}/{hc}.chr{chr}.bam",
        bqsr_report = rules.get_bqsr_report.output.bqsr_report,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        recal_bam = "data/recal_bams/{hc}.chr{chr}.recal.bam",
        recal_bai = "data/recal_bams/{hc}.chr{chr}.recal.bam.bai"
    shell: """
        gatk --java-options "-Xmx8G" ApplyBQSR \
        -I {input.dedup_bam_chunk} \
        -R {input.reference} \
        --bqsr-recal-file {input.bqsr_report} \
        -O {output.recal_bam}

        samtools index {output.recal_bam}
    """

rule haplotype_call:
    input:
        recal_bam = rules.apply_bqsr.output.recal_bam,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        gvcf = "results/call/vcfs/regions/{hc}/{hc}.chr{chr}.gvcf.vcf.gz"
    resources: mem = '20G'
    threads: 8
    shell: """
        gatk --java-options "-Xmx20G" HaplotypeCaller \
        -R {input.reference} \
        -I {input.recal_bam} \
        -O {output.gvcf} \
        -ERC GVCF
    """

rule genomics_db_import:
    input:
        gvcfs = expand("results/call/vcfs/regions/{hc}/{hc}.chr{chr}.gvcf.vcf.gz", hc = test_hc, allow_missing = True),
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        genomics_db = directory("results/call/tmp/chr{chr}.combined.db")
    resources:
        mem = '15G'
    threads: 4
    run:
        gvcf_files = ""
        for gvcf in input.gvcfs:
            gvcf_files = gvcf_files + "--variant " + gvcf + " "
        shell("""
        gatk --java-options "-Xmx55g -Xms2g" GenomicsDBImport \
        -R {input_ref} \
        -L {chromosome} \
        {gvcfs} \
        --tmp-dir "results/call/tmp/" \
        --genomicsdb-workspace-path {output}
        """.format(gvcfs = gvcf_files, input_ref = input.reference, chromosome = "chr" + wildcards.chr, output = output.genomics_db))

rule genotype_gvcf:
    input:
        gvcfs = rules.genomics_db_import.output.genomics_db,
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict
    output:
        called = "results/call/vcfs/regions/chr{chr}.vcf.gz",
        index = "results/call/vcfs/regions/chr{chr}.vcf.gz.tbi"
    resources:
        mem = '15G'
    threads: 4
    shell: """
        gatk --java-options "-Xmx55g" GenotypeGVCFs  \
        -R {input.reference} \
        -V gendb://{input.gvcfs} \
        -O {output.called}

        bcftools index -t {output.called}
    """

rule concat_hc_vcfs:
    input:
        vcfs = expand("results/call/vcfs/regions/chr{chr}.vcf.gz", chr = chromosomes),
        tbis = expand("results/call/vcfs/regions/chr{chr}.vcf.gz.tbi", chr = chromosomes)
    output:
        vcf = "results/call/vcfs/{hc}.vcf.gz" # Need to expand this if we're calling at different sites
    threads: 4
    resources: mem = '10G'
    shell: """
        bcftools concat --ligate -Oz -o {output.vcf} {input.vcfs}
    """

### variant calling END ###