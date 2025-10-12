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

samples_lc = read_tsv_as_lst(config['samples_lc'])
panels = config["panels"]

samples_chip = read_tsv_as_lst(config['samples_chip'])
seq_to_extract = [sample for sample in samples_lc if sample in samples_chip]

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

rule haplotype_call:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        ref_vcf = f"data/ref_panel/{hc_panel}/{hc_panel}.chr{{chr}}.vcf.gz" # Move this to config so that we can test sites at different ref panels
    output:
        vcf = f"results/call/vcfs/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz",
        empty_vcf1 = temp("results/call/tmp/ref/empty1_{type}_chr{chr}.vcf.gz"),
        empty_vcf2 = temp("results/call/tmp/ref/empty2_{type}_chr{chr}.vcf.gz")
    resources: mem = '100G'
    params:
        padding = 300
    threads: 16
    shell: """
        mkdir -p results/call/tmp/ref/
        mkdir -p results/call/vcfs/{hc_panel}/
        file=$(head -n 1 {input.bamlist})

        bcftools view -G {input.ref_vcf} | bcftools view -v {wildcards.type} -Oz -o {output.empty_vcf1}
        gatk IndexFeatureFile -I {output.empty_vcf1}

        gatk UpdateVCFSequenceDictionary \
        -V {output.empty_vcf1} \
        --source-dictionary $file \
        --output {output.empty_vcf2} \
        --replace true
        
        if [[ {wildcards.type} == "indels" ]]
        then
            gatk --java-options "-Xmx20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {output.empty_vcf2} \
            --alleles {output.empty_vcf2} \
            --native-pair-hmm-threads 16 \
            --assembly-region-padding {params.padding} \
            --output-mode EMIT_VARIANTS_ONLY
        else
            gatk --java-options "-Xmx20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {output.empty_vcf2} \
            --alleles {output.empty_vcf2} \
            --native-pair-hmm-threads 16 \
            --output-mode EMIT_VARIANTS_ONLY
        fi

        rm "{output.empty_vcf1}.tbi" "{output.empty_vcf2}.tbi"
    """

rule get_vqsr_report:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.haplotype_call.output.vcf
    output:
        tranch = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.tranches",
        recal = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.recal"
    params:
        hapmap = "data/GATK_resource_bundle/hapmap_3.3.hg38.vcf.gz",
        omni = "data/GATK_resource_bundle/1000G_omni2.5.hg38.vcf.gz",
        oneKG_snps = "data/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        GRCh38_indels = "data/GATK_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz",
        oneKG_indels = "data/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "data/GATK_resource_bundle/dbsnp_146.hg38.vcf.gz"
    resources:
        mem = '30G'
    shell: """
        mkdir -p results/call/VQSR/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.oneKG_snps} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=7 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP -O {output.recal} --tranches-file {output.tranch}
        else
            gatk --java-options "-Xms4G -Xmx4G" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.oneKG_indels} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL -O {output.recal} --tranches-file {output.tranch}
        fi
    """

rule apply_vqsr:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.haplotype_call.output.vcf,
        tranch = rules.get_vqsr_report.output.tranch,
        recal = rules.get_vqsr_report.output.recal
    output:
        recal_vcf = f"results/call/recal_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    resources:
        mem = '10G'
    shell: """
        mkdir -p results/call/recal_vcf/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.recal_vcf}
        else
            gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.recal_vcf}
        fi
    """

### variant calling END ###

rule get_chip_vcf:
    input:
        chip_qced = "results/chip/vcf/chip_qced.vcf.gz"
    output:
        chip_vcf = temp("results/chip/tmp/{id}/{id}.vcf.gz")
    params:
        seq_name = lambda wildcards: wildcards.id,
        chip_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Chip_Name'].values[0]
    resources:
        mem_mb = 30000
    shell: """
        mkdir -p results/chip/tmp/{wildcards.id}/

        if (bcftools query -l {input.chip_qced} | grep -q {params.chip_name}); then
            bcftools view -s {params.chip_name} -Oz -o {output.chip_vcf} {input.chip_qced}
        else
            touch {output.chip_vcf}
        fi
    """

rule get_imputation_vcf:
    input:
        imputation_result = "results/imputation/vcfs/{panel}/quilt.chr{chr}.vcf.gz"
    output:
        imputation_vcf = temp("results/imputation/tmp/{id}/{panel}_chr{chr}.vcf.gz")
    params:
        seq_name = lambda wildcards: wildcards.id,
        sample_name = lambda wildcards: sample_linker[sample_linker['Seq_Name'] == wildcards.id]['Sample_Name'].values[0]
    resources:
        mem_mb = 30000
    shell: """
        if (bcftools query -l {input.imputation_result} | grep -q {params.sample_name}); then
            bcftools view -s {params.sample_name} -Oz -o {output.imputation_vcf} {input.imputation_result}
        else
            touch {output.imputation_vcf}
        fi
    """

#vcf_dict = {}
#for panel in panels:
#    for id in seq_names:
#        vcf_dict[panel + id] = ["results/imputation/tmp/"+id+"/"+panel+"_chr"+str(chr)+".vcf.gz" for chr in chromosomes]

rule calculate_imputation_accuracy:
    input:
        imputation_vcf = expand("results/imputation/tmp/{id}/{panel}_chr{chr}.vcf.gz", chr = chromosome, allow_missing=True),
        chip_vcf = "results/chip/tmp/{id}/{id}.vcf.gz",
        afs = expand("data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt", chr = chromosome)
    output:
        r2 = "results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv"
    resources:
        mem_mb = 300000
    threads: 4
    run:
        chromosomes = chromosome
        vcfs = input.imputation_vcf
        mafs = input.afs
        chip = input.chip_vcf
        MAF_ary = np.array([0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.95, 1])

        check_lst = vcfs + [chip]
        for i in check_lst: # If there is an empty file, report all imputation accuracy to be -1
            if os.path.getsize(i) == 0:
                r2 = -1*np.ones((2, np.size(MAF_ary) - 1))
                r2 = pd.DataFrame(r2.T, columns = ['Imputation Accuracy','Bin Count'], index = MAF_ary[1:])
                r2.index.name = 'MAF'
                r2.to_csv(output.r2, sep=',', mode='a')
                sys.exit(0)

        vcf = lcwgsus.multi_parse_vcf(chromosomes, vcfs)
        af = lcwgsus.multi_read_af(chromosomes, mafs)
        chip = lcwgsus.read_vcf(chip)
        chip = lcwgsus.drop_cols(chip, drop_lst = ['id', 'qual', 'filter','info','format'])
        r2 = lcwgsus.calculate_imputation_accuracy(vcf, chip, af)
        r2.to_csv(output.r2, sep=',', mode='a')

rule plot_imputation_accuracy:
    input:
        r2 = expand("results/imputation/imputation_accuracy/{id}/{panel}_imputation_accuracy.csv", id = seq_to_extract, panel = panels)
    output:
        graph = "graphs/imputation_vs_chip.png"
    resources:
        mem_mb = 30000
    params:
        samples = seq_to_extract,
        panels = panels,
        linker = config['sample_linker']
    run:
        samples = params.samples
        panels = params.panels
        linker = pd.read_table(params.linker, sep = ',')
        fv = linker[(linker['Seq_Name'].isin(samples)) & (~linker['Sample_Name'].str.contains('mini'))]
        mini = linker[(linker['Seq_Name'].isin(samples)) & (linker['Sample_Name'].str.contains('mini'))]
        samples_fv = fv['Seq_Name'].to_list()
        samples_mini = mini['Seq_Name'].to_list()
        r2_fv, _ = lcwgsus.read_r2(panels, samples_fv)
        r2_fv = lcwgsus.aggregate_r2(r2_fv)
        r2_mini, _ = lcwgsus.read_r2(panels, samples_mini)
        r2_mini = lcwgsus.aggregate_r2(r2_mini)

        plt.figure(figsize = (10,6))
        for df in r2_fv:
            panel = df['panel'].values[0]
            plt.plot(np.arange(1, df.shape[0]+1), df['corr'], label = panel)
        for df in r2_mini:
            panel = df['panel'].values[0]
            plt.plot(np.arange(1, df.shape[0]+1), df['corr'], label = panel + '_mini', ls ='--')
        plt.xticks(np.arange(1, r2_fv[0].shape[0]+1), r2_fv[0]['AF'], rotation = 45)
        plt.xlabel('Allele frequencies (%)')
        plt.legend()
        plt.text(x = -1.5, y = 1.04, s = 'Aggregated imputation accuracy ($r^2$)')
        plt.grid(alpha = 0.5)
        plt.savefig(output.graph, bbox_inches = "tight", dpi=300)

rule calculate_imputation_accuracy_all:
    input:
        quilt_vcf = "results/imputation/vcfs/oneKG/quilt.chr{chr}.vcf.gz",
        chip_vcf = "results/chip/vcf/chip_by_chr/chip.chr{chr}.vcf.gz",
        af = "data/gnomAD_MAFs/afr/gnomAD_MAF_afr_chr{chr}.txt"
    output:
        h_report = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.tsv",
        h_impacc = "results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv",
        v_report = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.tsv",
        v_impacc = "results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv"
    resources:
        mem = '60G'
    threads: 8
    params:
        linker = config['sample_linker']
    run:
        mini = False
        common_cols = ['chr', 'pos', 'ref', 'alt']
        lc_sample_prefix = 'GM'
        chip_sample_prefix = 'GAM'
        seq_sample_prefix = 'IDT'

        quilt_vcf = input.quilt_vcf
        chip_vcf = input.chip_vcf
        af = input.af

        chip, lc, af = lcwgsus.imputation_calculation_preprocess(chip_vcf, quilt_vcf, af_txt)
        
        h_report = lcwgsus.calculate_h_imputation_accuracy(chip, lc, af, 
                                                   save_file = True, 
                                                   outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                                   save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.tsv')
        h_report = h_report.drop(columns = common_cols)

        h_impacc = lcwgsus.generate_h_impacc(h_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_variant/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.h.impacc.tsv')
                                           
        v_report = lcwgsus.calculate_v_imputation_accuracy(chip, lc, af, 
                                           save_file = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.tsv')

        v_impacc = lcwgsus.generate_v_impacc(v_report, 
                                           save_impacc = True, 
                                           outdir = 'results/imputation_metrics/lc_chip/all_samples/by_sample/', 
                                           save_name = 'lc.chip.typed.chr' + str(chromosome) +'.v.impacc.tsv')
        # Ignore the _BC fields in vertical reports as they are not reliable

rule plot_imputation_accuracy:
    input:
        h_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_variant/lc.chip.typed.chr{chr}.h.impacc.tsv", chr = chromosome),
        v_impaccs = expand("results/imputation_metrics/lc_chip/all_samples/by_sample/lc.chip.typed.chr{chr}.v.impacc.tsv", chr = chromosome)
    output:
        r2NRC_h = "graphs/imputation/lc_chip/all_samples/by_variant/r2_NRC.png",
        ccd_h = "graphs/imputation/lc_chip/all_samples/by_variant/ccd_by_genotype.png",
        r2NRC_v = "graphs/imputation/lc_chip/all_samples/by_sample/r2_NRC.png",
        ccd_v = "graphs/imputation/lc_chip/all_samples/by_sample/ccd_by_genotype.png"
    resources:
        mem = '30G'
    run:
        h_lst = input.h_impaccs
        v_lst = input.v_impaccs
        outdir_h = "graphs/imputation/lc_chip/all_samples/by_variant/"
        outdir_v = "graphs/imputation/lc_chip/all_samples/by_sample/"

        h_dfs = [pd.read_csv(i, sep = '\t') for i in h_lst]
        h = lcwgsus.average_impacc_by_chr(h_dfs)
        v_dfs = [pd.read_csv(i, sep = '\t') for i in v_lst]
        v = lcwgsus.average_impacc_by_chr(v_dfs)

        dfs = [h[['AF', 'r2', 'r2_AC']], h[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by variant', save_fig = True, outdir = outdir_h, save_name = "r2_NRC.png")

        dfs = [h[['AF', 'ccd_homref', 'ccd_homref_AC']], h[['AF', 'ccd_het', 'ccd_het_AC']], h[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_h, save_name = "ccd_by_genotype.png")
    
        dfs = [v[['AF', 'r2', 'r2_AC']], v[['AF', 'NRC', 'NRC_AC']]]
        lcwgsus.plot_imputation_accuracy(dfs, title = 'r2 and NRC by sample', save_fig = True, outdir = outdir_v, save_name = "r2_NRC.png")

        dfs = [v[['AF', 'ccd_homref', 'ccd_homref_AC']], v[['AF', 'ccd_het', 'ccd_het_AC']], v[['AF', 'ccd_homalt', 'ccd_homalt_AC']]]
        plot_imputation_accuracy(dfs, title = 'Comparing different genotypes', save_fig = True, outdir = outdir_v, save_name = "ccd_by_genotype.png")


rule hla_clean_bam_alt:
    input:
        bam = rules.hla_imputation_preprocess_alt.output.chr
    output:
        bam = "data/hla_bams_alt/{id}.chr6.bam",
        bai = "data/hla_bams_alt/{id}.chr6.bam.bai",
        sam = temp("data/hla_bams_alt/{id}.chr6.sam"),
        tmp1 = temp("data/hla_bams_alt/{id}.tmp1.bam"),
        metric = temp("data/hla_bams_alt/{id}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/hla_bams_alt/tmp/{id}/",
        sample = "{id}"
    shell: """
        mkdir -p {params.tmpdir}

        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {input.bam}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        samtools index {output.bam}

        samtools view -h {output.bam} chr6:26000000-34000000 > {output.sam}
        samtools view {output.bam} | \
        awk -F '\t' '($3~/HLA/){{print}}' >> {output.sam}
        samtools view -bS {output.sam} > {output.tmp1}

        samtools view -H {output.bam} > {output.sam}
        samtools view {output.tmp1} | \
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1 ~ /^@/) {{
                print $0
            }} else {{
                new_qual = ""
                qual = $11
                for (i = 1; i <= length(qual); i++) {{
                    q = substr(qual, i, 1)
                    if (q ~ /[@ABCDEF]/) {{
                        q = "?"
                    }}
                    new_qual = new_qual q
                }}
                $11 = new_qual
                print $0
            }}
        }}' >> {output.sam}
        samtools view -bS {output.sam} > {output.bam}

        samtools index {output.bam}
    """

rule alignment:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
    output:
        bam = temp("data/bams/tmp/{id}.bam")
    resources:
        mem = '10G'
    threads: 8
    params:
        reheader = reheader
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
        if [[ -d "data/bam_headers" && {params.reheader} == "True" ]]
        then
            mkdir -p data/tmp
            grep "^@RG" "data/bam_headers/{wildcards.id}.header.txt" | tr '\t' '\n' > "data/tmp/{wildcards.id}.rgs.txt"
            picard AddOrReplaceReadGroups \
            I={output.bam} \
            O="data/tmp/{wildcards.id}.bam" \
            RGLB=$(grep "^LB:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/LB://g') \
            RGPL=$(grep "^PL:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/PL://g') \
            RGPU=$(grep "^PU:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/PU://g') \
            RGSM=$(grep "^SM:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/SM://g') \
            RGID=$(grep "^ID:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/ID://g') \
            RGCN=$(grep "^CN:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/CN://g') \
            RGDT=$(grep "^DT:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/DT://g')
            rm {output.bam} "data/tmp/{wildcards.id}.rgs.txt"
            mv "data/tmp/{wildcards.id}.bam" {output.bam}
        fi
    """

# if clean_fastq:
#     ruleorder: alignment_alt > alignment
# else:
#     ruleorder: alignment > alignment_alt

rule alignment_alt:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
    output:
        bam = temp("data/bams/tmp/{id}.bam")
    resources:
        mem = '10G'
    threads: 8
    params:
        reheader = reheader
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.bam}
        if [[ -d "data/bam_headers" && {params.reheader} == "True" ]]
        then
            mkdir -p data/tmp
            grep "^@RG" "data/bam_headers/{wildcards.id}.header.txt" | tr '\t' '\n' > "data/tmp/{wildcards.id}.rgs.txt"
            picard AddOrReplaceReadGroups \
            I={output.bam} \
            O="data/tmp/{wildcards.id}.bam" \
            RGLB=$(grep "^LB:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/LB://g') \
            RGPL=$(grep "^PL:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/PL://g') \
            RGPU=$(grep "^PU:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/PU://g') \
            RGSM=$(grep "^SM:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/SM://g') \
            RGID=$(grep "^ID:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/ID://g') \
            RGCN=$(grep "^CN:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/CN://g') \
            RGDT=$(grep "^DT:" "data/tmp/{wildcards.id}.rgs.txt" | sed 's/DT://g')
            rm {output.bam} "data/tmp/{wildcards.id}.rgs.txt"
            mv "data/tmp/{wildcards.id}.bam" {output.bam}
        fi
    """

rule fixmate:
    input:
        bam = rules.alignment_alt.output.bam if clean_fastq else rules.alignment.output.bam
    output:
        fixmate = temp("data/bams/tmp/{id}.fixmate.bam")
    resources: mem = '10G'
    shell: """
        samtools sort -n {input.bam} | samtools fixmate -m - {output.fixmate}
    """

rule sort:
    input:
        fixmate = rules.fixmate.output.fixmate
    output:
        sorted = temp("data/bams/tmp/{id}.sorted.bam")
    resources: mem = '10G'
    shell: """
        samtools sort -o {output.sorted} {input.fixmate}
    """

rule index:
	input:
		sorted = rules.sort.output.sorted
	output:
		bai = temp("data/bams/tmp/{id}.sorted.bam.bai")
	resources: mem = '10G'
	shell: """
		samtools index {input.sorted}
	"""

rule rename_alignment_files:
	input:
		bam = rules.index.input.sorted,
		bai = rules.index.output.bai
	output:
		bam = "data/bams/{id}.bam",
		bai = "data/bams/{id}.bam.bai"
	shell: """
		mv {input.bam} {output.bam}
		mv {input.bai} {output.bai}
	"""

ule hla_clean_bam:
    input:
        bam = rules.hla_imputation_preprocess.output.chr
    output:
        bam = "data/hla_bams/{id}.chr6.bam",
        bai = "data/hla_bams/{id}.chr6.bam.bai",
        sam = temp("data/hla_bams/{id}.chr6.sam"),
        tmp1 = temp("data/hla/bams/{id}.tmp1.bam"),
        metric = temp("data/hla/bams/{id}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/hla/bams/tmp/{id}/",
        sample = "{id}"
    shell: """
        mkdir -p {params.tmpdir}

        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {input.bam}

        picard MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.bam} {output.tmp1}

        samtools index {output.bam}

# Filter chr6 and HLA contigs as well as recode QUAL strs

        samtools view -h {output.bam} chr6:26000000-34000000 > {output.sam}
        samtools view {output.bam} | \
        awk -F '\t' '($3~/HLA/){{print}}' >> {output.sam}
        samtools view -bS {output.sam} > {output.tmp1}

        samtools view -H {output.bam} > {output.sam}
        samtools view {output.tmp1} | \
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1 ~ /^@/) {{
                print $0
            }} else {{
                new_qual = ""
                qual = $11
                for (i = 1; i <= length(qual); i++) {{
                    q = substr(qual, i, 1)
                    if (q ~ /[@ABCDEF]/) {{
                        q = "?"
                    }}
                    new_qual = new_qual q
                }}
                $11 = new_qual
                print $0
            }}
        }}' >> {output.sam}
        samtools view -bS {output.sam} > {output.bam}

        samtools index {output.bam}
    """

rule hla_clean_bam_alt:
    input:
        bam = rules.hla_imputation_preprocess_alt.output.chr
    output:
        bam = "data/hla_bams_alt/{id}.chr6.bam",
        bai = "data/hla_bams_alt/{id}.chr6.bam.bai",
        tmp1 = temp("data/hla_bams_alt/{id}.tmp1.bam"),
        metric = temp("data/hla_bams_alt/{id}.metrics.txt")
    threads: 8
    resources:
        mem = '50G'
    params:
        tmpdir = "data/hla_bams_alt/tmp/{id}/",
        sample = "{id}",
        picard = tools['picard_plus']
    shell: """
        mkdir -p {params.tmpdir}

        picard FixMateInformation -I {input.bam}

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.tmp1} {input.bam}

        {params.picard} MarkDuplicates \
        -I {output.tmp1} \
        -O {output.bam} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -T {params.tmpdir} -o {output.tmp1} {output.bam}

        samtools index {output.tmp1}

        samtools view -h {output.tmp1} -o {output.bam} chr6
        samtools index {output.bam}
    """


rule index:
    input:
        reference = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        indexed = "data/references/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb"
    resources:
        mem = '50G'
    threads: 8
    shell: """
        bwa index {input.reference}
    """

def get_num_reads(wildcards):
    total = 3200000000
    cov = 0.6
    num = round(total*cov/(2*int(wildcards.rl.split('-')[0])))
    return num

rule simulate_reads:
    input:
        fasta = "data/lr_fasta/HG02886.{hap}.fa"
    output:
        fastq1 = temp("data/sr_simulations/{rl}/tmp.{hap}.{rl}.bwa.read1.fastq.gz"),
        fastq2 = temp("data/sr_simulations/{rl}/tmp.{hap}.{rl}.bwa.read2.fastq.gz")
    resources:
        mem = '30G'
    threads: 4
    params:
        num_reads = get_num_reads,
        read_length = get_read_length,
        mean_length = get_mean_length,
        sd_length = get_mean_std,
        error1 = get_error1,
        error2 = get_error2,
        outdir = "data/sr_simulations/{rl}/",
        output_prefix = "data/sr_simulations/{rl}/tmp.{hap}.{rl}"
    shell: """
        mkdir -p {params.outdir}
        
        dwgsim -N {params.num_reads} \
        -1 {params.read_length} \
        -2 {params.read_length} \
        -e {params.error1} -E {params.error2} \
        -d {params.mean_length} \
        -s {params.sd_length} \
        -y 0 \
        {input.fasta} {params.output_prefix}
    """

# ";NA;" is a common feature for a match in both ref panels
rule merge_malariaGen_v3_with_oneKG:
    input:
        mg = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz",
        oneKG = "data/ref_panel/oneKG/oneKG.chr{chr}.vcf.gz"
    output:
        vcf = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz",
        tbi = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz.tbi"
    resources: mem = '70G'
    threads: 4
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b38/

        tabix -f {input.mg}
        
        bcftools merge {input.mg} {input.oneKG} -Ou | \
        bcftools view -i 'ID~";NA;"' -Oz -o {output.vcf}

        tabix -f {output.vcf} 
    """

rule vanilla_alignment:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta"
    output:
        bam = "data/bams/{id}.bam",
        bai = "data/bams/{id}.bam.bai",
        tmp1 = temp("data/bams/{id}.tmp1.bam"),
        metric = temp("data/bams/{id}.metrics.txt")
    resources:
        mem = '30G'
    params: 
        sample = "{id}",
        picard = tools["picard_plus"]
    threads: 6
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.tmp1}
        
        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}

        samtools index {output.bam}

        picard AddOrReplaceReadGroups \
        -VERBOSITY ERROR \
        -I {output.bam} \
        -O {output.tmp1} \
        -RGLB OGC \
        -RGPL Illumina \
        -RGPU unknown \
        -RGSM {params.sample}

        {params.picard} FixMateInformation -I {output.tmp1}

        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}

        {params.picard} MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

rule merge_malariaGen_v3_with_oneKG:
    input:
        mg = "data/ref_panel/malariaGen_v3_b38_alone/malariaGen_v3_b38_alone.chr{chr}.vcf.gz",
        oneKG = "data/ref_panel/oneKG/oneKG.chr{chr}.vcf.gz"
    output:
        vcf = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz",
        tbi = "data/ref_panel/malariaGen_v3_b38/malariaGen_v3_b38.chr{chr}.vcf.gz.tbi"
    resources: mem = '70G'
    threads: 4
    shell: """
        mkdir -p data/ref_panel/malariaGen_v3_b38/

        tabix -f {input.mg}
        
        bcftools merge -0 {input.mg} {input.oneKG} -Oz -o {output.vcf}

        tabix -f {output.vcf} 
    """

rule subset_lc_samples:
    input:
        quilt_vcf = '{imp_dir}vcf/all_samples/lc_vcf/lc.chr{chr}.vcf.gz',
        chip_vcf = '{imp_dir}vcf/all_samples/hc_vcf/hc.chr{chr}.vcf.gz'
    output:
        ss_vcf = temp('{imp_dir}vcf/all_samples/lc_vcf/lc.subset.chr{chr}.vcf.gz'),
        tmp_names = temp('{imp_dir}vcf/all_samples/samples.chr{chr}.tsv')
    resources:
        mem = '10G'
    run: 
        hc_names = lcwgsus.bcftools_get_samples(input.chip_vcf)
        lc_names = lcwgsus.bcftools_get_samples(input.quilt_vcf)

        if lc_names[0].startswith('IDT'):
            sl = sample_linker.sort_values(by = 'Seq_Name')
            sl = sl[sl['Seq_Name'].isin(lc_names)]
            sl = sl[~sl['Sample_Name'].str.contains('mini')]
            rename_map = {k:v for k,v in zip(sl['Chip_Name'], sl['Seq_Name'])}
            samples = []
            for n in hc_names:
                samples.append(rename_map[n])
            lcwgsus.save_lst(output.tmp_names, samples)

        elif lc_names[0].startswith('GM'):
            rename_map = lcwgsus.generate_rename_map()
            lc = wildcards.imp_dir.split('/')[-2].split('_')[0]
            samples = lcwgsus.find_matching_samples(hc_names, rename_map, lc = lc)
            lcwgsus.save_lst(output.tmp_names, samples)
        else:
            rename_map = lcwgsus.generate_rename_map()
            lc = wildcards.imp_dir.split('/')[-2].split('_')[0]
            samples = lcwgsus.find_matching_samples(hc_names, rename_map, lc = lc)
            lcwgsus.save_lst(output.tmp_names, samples)

        shell("bcftools view -S {output.tmp_names} -Oz -o {output.ss_vcf} {input.quilt_vcf}")

REGIONS={}
for chr in chromosome:
    start=[10000001, 15000001]
    end=[  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

file="results/imputation/regions.json"
if os.path.exists(file):
    print("Replacing regions to impute with derived file")
    with open(file) as json_file:
        REGIONS = json.load(json_file)

regions_to_prep = []
vcfs_to_impute = []
vcfs_to_concat = {}
for chr in chromosome:
    start = REGIONS[str(chr)]["start"]
    end = REGIONS[str(chr)]["end"]
    per_chr_vcfs = []
    for i in range(0, start.__len__()):
        regionStart = start[i]
        regionEnd = end[i]
        file = "results/imputation/refs/" + PANEL_NAME + "/RData/ref_package.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".RData"
        regions_to_prep.append(file)
        file = "results/imputation/vcfs/" + PANEL_NAME + "/regions/quilt.chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
        vcfs_to_impute.append(file)
        per_chr_vcfs.append(file)
    vcfs_to_concat[str(chr)] = per_chr_vcfs
    return regions_to_prep, vcfs_to_impute, vcfs_to_concat

rule fix_bam:
    input:
        bam = "data/bams/{id}.bam"
    output:
        tmp1 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.tmp1.bam"),
        tmp2 = temp("data/chunk_bams/tmp/tmp/{id}/{id}.tmp2.bam")
    threads: 4
    resources: mem = '10G'
    params:
        verbosity = "ERROR",
        sample = "{id}".split("_")[0]
    shell: """
        mkdir -p data/chunk_bams/tmp/tmp/{wildcards.id}/
        samtools sort -o {output.tmp1} {input.bam}

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
    """

REGIONS = {}
for chr in chromosome:
    start = [10000001, 15000001]
    end = [  15000000, 20000000]
    REGIONS[str(chr)] = {"start":start, "end":end}

chunk_region = "results/imputation/regions.json"
if os.path.exists(chunk_region):
    with open(chunk_region) as json_file:
        REGIONS = json.load(json_file)

region_vcfs = []
for chr in chromosome:
    start = REGIONS[str(chr)]["start"]
    end = REGIONS[str(chr)]["end"]
    for t in variant_types:
        for i in range(0, start.__len__()):
            regionStart = start[i]
            regionEnd = end[i]
            file = "results/call/recal_vcf/" + hc_panel + "/regions/" + hc_panel + "." + t +  ".chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            region_vcfs.append(file)

rule variant_calling_all:
    input:
        fai = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta.fai" if concatenate else "data/references/GRCh38.fa.fai",
        dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict" if concatenate else "data/references/GRCh38.dict",
        bqsr_known_sites = [file + ".tbi" for file in config["bqsr_known_sites"]],
        gatk_to_index = [file + ".tbi" for file in config["gatk_to_index"]],
        #bqsr_reports = expand("results/call/BQSR/BQSR_reports/{hc}.BQSR.report", hc = samples_hc),
        #recal_bams = expand("data/recal_bams/{hc}.recal.bam", hc = samples_hc),
        #recal_bais = expand("data/recal_bams/{hc}.recal.bam.bai", hc = samples_hc),
        #bamlist = "results/call/bam.list",
        regions = [region_vcfs],
        merge_vcf = expand(f"results/call/merge_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz", type = variant_types, chr = chromosome),
        tranches = expand(f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.tranch", type = variant_types, chr = chromosome),
        recals = expand(f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.recal", type = variant_types, chr = chromosome),
        recal_vcf = expand(f"results/call/recal_vcf/{hc_panel}/{{type}}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz", type = variant_types, chr = chromosome)

rule GATK_chunk_reference:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        ref_vcf = f"data/ref_panel/{hc_panel}/{hc_panel}.chr{{chr}}.vcf.gz" # Move this to config so that we can test sites at different ref panels
    output:
        empty_vcf1 = temp("results/call/tmp/ref/empty1_{type}_{regionStart}.{regionEnd}_chr{chr}.vcf.gz"),
        empty_vcf2 = temp("results/call/tmp/ref/empty2_{type}_{regionStart}.{regionEnd}_chr{chr}.vcf.gz")
    resources: mem = '10G'
    shell: """
        mkdir -p results/call/tmp/ref/
        mkdir -p results/call/vcfs/{hc_panel}/
        file=$(head -n 1 {input.bamlist})

        bcftools view -G -v {wildcards.type} -r chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} -Oz -o {output.empty_vcf1} {input.ref_vcf}
        
        gatk IndexFeatureFile -I {output.empty_vcf1}

        gatk UpdateVCFSequenceDictionary \
        -V {output.empty_vcf1} \
        --source-dictionary $file \
        --output {output.empty_vcf2} \
        --replace true
    """

rule haplotype_call:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        bamlist = rules.prepare_hc_bamlist.output.bamlist,
        empty_vcf2 = rules.GATK_chunk_reference.output.empty_vcf2
    output:
        vcf = f"results/call/vcfs/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz"
    resources: mem = '20G'
    params:
        padding = 300
    threads: 4
    shell: """
        if [[ {wildcards.type} == "indels" ]]
        then
            gatk --java-options "-Xmx20G -Xms20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {input.empty_vcf2} \
            --alleles {input.empty_vcf2} \
            --native-pair-hmm-threads 4 \
            --assembly-region-padding {params.padding} \
            --output-mode EMIT_VARIANTS_ONLY
        else
            gatk --java-options "-Xmx20G -Xms20G" HaplotypeCaller \
            -R {input.reference} \
            -I {input.bamlist} \
            -O {output.vcf} \
            -L {input.empty_vcf2} \
            --alleles {input.empty_vcf2} \
            --native-pair-hmm-threads 4 \
            --output-mode EMIT_VARIANTS_ONLY
        fi

        rm "{input.empty_vcf2}.tbi"
    """

REGIONS = {}
for chr in chromosome:
    start = [10000001, 15000001]
    end = [  15000000, 20000000]
    REGIONS[str(chr)]={"start":start, "end":end}

chunk_regions = "results/imputation/regions.json"
if os.path.exists(chunk_regions):
    with open(chunk_regions) as json_file:
        REGIONS = json.load(json_file)

vcfs_to_concat = {}
for chr in chromosome:
    start = REGIONS[str(chr)]["start"]
    end = REGIONS[str(chr)]["end"]
    vcfs_to_concat[str(chr)] = {}
    for t in variant_types:
        file_ary = []
        for i in range(0, start.__len__()):
            regionStart = start[i]
            regionEnd = end[i]
            file = "results/call/vcfs/" + hc_panel + "/" + hc_panel + "." + t +  ".chr" + str(chr) + "." + str(regionStart) + "." + str(regionEnd) + ".vcf.gz"
            file_ary.append(file)
        vcfs_to_concat[str(chr)][t] = file_ary

rule concat_hc_vcfs:
    input:
        vcfs = lambda wildcards: vcfs_to_concat[str(wildcards.chr)][str(wildcards.type)]
    output:
        vcf = f"results/call/merge_vcf/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz"
    threads: 4
    resources: mem = '20G'
    shell: """
        mkdir -p results/call/merge_vcf/{hc_panel}/

        bcftools concat -Oz -o {output.vcf} -a -d {wildcards.type} {input.vcfs}
        tabix {output.vcf}
    """

rule get_vqsr_report:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.concat_hc_vcfs.output.vcf
    output:
        tranch = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.tranch",
        recal = f"results/call/VQSR/{hc_panel}/{hc_panel}.{{type}}.chr{{chr}}.recal"
    params:
        hapmap = "data/GATK_resource_bundle/hapmap_3.3.hg38.vcf.gz",
        omni = "data/GATK_resource_bundle/1000G_omni2.5.hg38.vcf.gz",
        oneKG_snps = "data/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        GRCh38_indels = "data/GATK_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz",
        oneKG_indels = "data/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "data/GATK_resource_bundle/dbsnp_146.hg38.vcf.gz"
    resources:
        mem = '20G'
    shell: """
        mkdir -p results/call/VQSR/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.oneKG_snps} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=7 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP -O {output.recal} --tranches-file {output.tranch}
        else
            gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \
            -tranche 99.0 \
            -R {input.reference} \
            -V {input.vcf} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.oneKG_indels} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode INDEL -O {output.recal} --tranches-file {output.tranch}
        fi
    """

rule apply_vqsr:
    input:
        reference = rules.GATK_prepare_reference.input.reference,
        fai = rules.GATK_prepare_reference.output.fai,
        dict = rules.GATK_prepare_reference.output.dict,
        vcf = rules.concat_hc_vcfs.output.vcf,
        tranch = rules.get_vqsr_report.output.tranch,
        recal = rules.get_vqsr_report.output.recal
    output:
        recal_vcf = f"results/call/recal_vcf/{hc_panel}/{{type}}/{hc_panel}.{{type}}.chr{{chr}}.vcf.gz",
        tmp_vcf = temp(f"results/call/recal_vcf/{hc_panel}/{{type}}/tmp.{hc_panel}.{{type}}.chr{{chr}}.vcf.gz")
    resources:
        mem = '20G'
    params:
        rename_samples = config['hc_vcf_rename_samples']
    shell: """
        mkdir -p results/call/recal_vcf/{hc_panel}/

        if [[ {wildcards.type} == "snps" ]]
        then
            gatk --java-options "-Xmx20g -Xms20g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.recal_vcf}
        else
            gatk --java-options "-Xmx20g -Xms20g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.recal_vcf}
        fi
        
        bcftools norm -m-any -Oz -o {output.tmp_vcf} {output.recal_vcf}
        rm {output.recal_vcf}

        if ! [ -f {params.rename_samples} ]; then
            bcftools reheader -s {params.rename_samples} -o {output.recal_vcf} {output.tmp_vcf}
        fi

    """

rule rmdup:
    input:
        bam = "data/bams/{id}.bam"
    output:
        dedup_bam = "data/dedup_bams/{id}.bam"
    resources: mem = '10G'
    shell: """
        samtools rmdup {input.bam} {output.dedup_bam}
    """

rule index_dedup:
    input:
        dedup_bam = rules.rmdup.output.dedup_bam
    output:
        dedup_bai = "data/dedup_bams/{id}.bam.bai"
    resources: mem = '10G'
    shell: """
        samtools index {input.dedup_bam}
    """

rule convert_chr6_to_chip_form:
    input:
        vcf = "{three_stage_vcf_outdir}chr6.tmp.vcf.gz"
    output:
        vcf = "{three_stage_vcf_outdir}chr6.vcf.gz"
    resources:
        mem = '30G'
    threads: 4
    run:
        imp_vcf = input.vcf
        lc = lcwgsus.read_vcf(imp_vcf)
        metadata = lcwgsus.read_metadata(imp_vcf)

        lc = lc.apply(lcwgsus.convert_to_chip_format, axis = 1)

        lcwgsus.save_vcf(lc,
             metadata,
             prefix='chr',
             outdir=wildcards.three_stage_vcf_outdir,
             save_name="chr6.vcf.gz"
             )
        lcwgsus.rezip_vcf(output.vcf)

rule fastuniq: # Currently deprecated as we are basically gonna remove in markdup
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        fastq1 = temp("data/tmp/{id}_fast_uniq_1.fastq.gz"),
        fastq2 = temp("data/tmp/{id}_fast_uniq_2.fastq.gz"),
        filelist = temp("data/tmp/{id}.list"),
        fastq1_uncompress = temp("data/tmp/input_{id}_fast_uniq_1.fastq.gz"),
        fastq2_uncompress = temp("data/tmp/input_{id}_fast_uniq_2.fastq.gz"),
        fastq1_unzip = temp("data/tmp/{id}_fast_uniq_1.fastq"),
        fastq2_unzip = temp("data/tmp/{id}_fast_uniq_2.fastq")
    threads: 4
    resources:
        mem = '30G'
    shell: """
        pigz -p {threads} -dc {input.fastq1} > {output.fastq1_uncompress}
        pigz -p {threads} -dc {input.fastq2} > {output.fastq2_uncompress}
        echo {output.fastq1_uncompress} > {output.filelist}
        echo {output.fastq2_uncompress} >> {output.filelist}
        fastuniq -i {output.filelist} -o {output.fastq1_unzip} -p {output.fastq2_unzip}
        gzip -c {output.fastq1_unzip} > {output.fastq1}
        gzip -c {output.fastq2_unzip} > {output.fastq2}
    """

hla_ref_panel_start = config["hla_b38_start_extra"]
hla_ref_panel_end = config["hla_b38_end_extra"]

rule merge_1KG_GAMCC_hla_only:
    input:
        haps = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr6.hap", to_merge = to_merge),
        legends = expand("results/hla_ref_panel/oneKG_mGenv1/tmp/{to_merge}.chr6.legend", to_merge = to_merge),
        gen_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr6-final.b38.txt",
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

rule phase_1KG_alleles:
    input:
        samples_file = "results/hla/imputation/ref_panel/auxiliary_files/oneKG.samples",
        hlatypes_file = "results/hla/imputation/ref_panel/auxiliary_files/20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        phased_vcf_file = "/well/band/users/rbx225/recyclable_files/ref_panels/oneKG_30x/oneKG.chr6.vcf.gz"
    output:
        # html = "results/phasing/html/oneKG-{filter}-{gene}.html",
        phase_df = "results/phasing/phased_dfs/oneKG-{filter}-{gene}.tsv"
    resources:
        mem = '30G'
    threads: 4
    params:
        ipd_gen_file_dir = '/well/band/users/rbx225/recyclable_files/hla_reference_files/alignments/',
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv',
        reference_allele_file = '/well/band/users/rbx225/recyclable_files/hla/b38_reference_alleles.tsv'
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        hlatypes = pd.read_csv(input.hlatypes_file, sep = '\t')

        ref_samples = pd.read_csv(input.samples_file, sep = ' ')
        ref_samples_removed = ref_samples[~ref_samples['SAMPLE'].isin(hlatypes['Sample ID'].tolist())]
        samples_to_remove = ref_samples_removed['SAMPLE'].tolist()
        hlatypes = hlatypes[~hlatypes['Sample ID'].isin(samples_to_remove)].sort_values(by = 'Sample ID').reset_index(drop = True)

        reference_allele_ary = np.array(lcwgsus.read_tsv_as_lst(params.reference_allele_file))
        strict_snp_filter = True if wildcards.filter == 'strict' else False

        return_dict = phase_hla_on_haplotypes(gene = wildcards.gene, 
                                    ipd_gen_file_dir = params.ipd_gen_file_dir, 
                                    hla_gene_information = hla_gene_information,
                                    hlatypes = hlatypes,
                                    phased_vcf = input.phased_vcf_file, 
                                    reference_allele_ary = reference_allele_ary, 
                                    strict_snp_filter = strict_snp_filter,
                                    read_from_QUILT = False, 
                                    subset_vcf_samples = None,
                                    sample_linker = None)
        # hlatypes = return_dict['hlatypes']

        # individual = 'NA12878'
        # ix = hlatypes.index[hlatypes['Sample ID'] == individual][0]
        # display_indices = np.arange(10)

        # res = visualise_phase(wildcards.gene, ix, hlatypes, return_dict, both_het = True)
        # compare_phase(display_indices, res, save_html = True, save_name = output.html)
        df = return_dict['phase_df'][['Sample', 'allele1', 'allele2']]
        df.columns = ['Sample ID', f'HLA-{wildcards.gene} 1', f'HLA-{wildcards.gene} 2']
        df.to_csv(output.phase_df, sep = '\t', index = False, header = True)


rule calculate_phasing_concordance:
    input:
        unphased_vcf = "results/phasing/HLA_1KG_BEAGLE/unphased.1KG.chr6.vcf.gz",
        phased_vcf = "results/phasing/HLA_1KG_BEAGLE/phased.1KG.chr6.vcf.gz",
        extracted_vcf = "results/phasing/HLA_1KG_BEAGLE/tmp/phased.{gene}.1KG.chr6.vcf.gz",
        phase_df = "results/phasing/phased_dfs/oneKG_{vcf_version}-{filter}-{gene}.tsv"
    output:
        concordance_df = "results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}-{gene}.tsv"
    params:
        hla_gene_information_file = '/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv'
    resources: mem = '80G'
    threads: 8
    run: 
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        our_hla = pd.read_csv(input.phase_df, sep = '\t')
        
        gene = wildcards.gene
        start = hla_gene_information[hla_gene_information['Name'] == f'HLA-{gene}']['Start'].values[0]
        end = hla_gene_information[hla_gene_information['Name'] == f'HLA-{gene}']['End'].values[0]

        vcf = lcwgsus.read_vcf(input.extracted_vcf)
        vcf = vcf.drop(columns = ['chr', 'pos', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT']).reset_index(drop = True).T

        beagle_hla = our_hla.copy()
        beagle_hla.iloc[:, 1:] = 'N/A'

        for i,s in enumerate(our_hla['Sample ID']):
            vcf_row = vcf.loc[s, :]
            alleles_idx = vcf_row.index[vcf_row.str.contains('1')]
            if len(alleles_idx) != 0:
                for j in alleles_idx:
                    name = vcf.loc['ID',j]
                    genotype = vcf.loc[s,j]
                    twofield = name.split('*')[1]
                    if genotype == '1|0':
                        beagle_hla.loc[i, f'HLA-{gene} 1'] = twofield
                    elif genotype == '0|1':
                        beagle_hla.loc[i, f'HLA-{gene} 2'] = twofield
                    elif genotype == '1|1':
                        beagle_hla.loc[i, f'HLA-{gene} 1'] = twofield
                        beagle_hla.loc[i, f'HLA-{gene} 2'] = twofield
                    else:
                        pass
        del vcf
        new = read_vcf(start = max(25000000, start - 10000), end = min(34000000, end + 10000), phased_vcf = input.phased_vcf, hlatypes = beagle_hla)
        old = read_vcf(start = max(25000000, start - 10000), end = min(34000000, end + 10000), phased_vcf = input.unphased_vcf, hlatypes = beagle_hla)
        new = pd.merge(new, old[['snp', 'pos']], how = 'inner')
        old = old.reset_index(drop = True)
        new_archive = new.copy()
        old_archive = old.copy()

        for i,s in enumerate(beagle_hla['Sample ID'].values):
            print(f'phasing sample {i}: {s}')
            new = new_archive.copy()
            old = old_archive.copy()
            tmp = new[new['pos']<start]
            idx_lst = tmp.index[tmp[s].isin(['0|1', '1|0'])]
            if len(idx_lst) != 0:
                idx = [idx_lst[-1]]
            else:
                tmp = new[new['pos']>end]
                idx_lst = tmp.index[tmp[s].isin(['0|1', '1|0'])]
                if len(idx_lst) != 0:
                    idx = [idx_lst[-1]]
                else:
                    idx = []

            multiplier = 2
            while (len(idx) == 0) or (multiplier > 5):
                tmp_new_up = read_vcf(start = max(25000000, start - multiplier*10000), end = max(25000000, start - (multiplier - 1)*10000), phased_vcf = input.phased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                tmp_old_up = read_vcf(start = max(25000000, start - multiplier*10000), end = max(25000000, start - (multiplier - 1)*10000), phased_vcf = input.unphased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                tmp_new_up = pd.merge(tmp_new_up, tmp_old_up[['snp', 'pos']], how = 'inner')
                tmp_old_up = tmp_old_up.reset_index(drop = True)
                
                idx_lst = tmp_new_up.index[tmp_new_up[s].isin(['0|1', '1|0'])]
                if len(idx_lst) != 0:
                    idx.append(idx_lst[-1])
                    new = tmp_new_up
                    old = tmp_old_up
                else:
                    tmp_new_down = read_vcf(start = min(34000000, end + (multiplier - 1)*10000), end = min(34000000, end + multiplier*10000), phased_vcf = input.phased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                    tmp_old_down = read_vcf(start = min(34000000, end + (multiplier - 1)*10000), end = min(34000000, end + multiplier*10000), phased_vcf = input.unphased_vcf, hlatypes = beagle_hla, subset_vcf_samples = s)
                    tmp_new_down = pd.merge(tmp_new_down, tmp_old_down[['snp', 'pos']], how = 'inner')
                    tmp_old_down = tmp_old_down.reset_index(drop = True)

                    idx_lst = tmp_new_down.index[tmp_new_down[s].isin(['0|1', '1|0'])]
                    if len(idx_lst) != 0:
                        idx.append(idx_lst[-1])
                        new = tmp_new_down
                        old = tmp_old_down

                multiplier += 1

            if len(idx) != 0:
                idx = idx[0]
                if new.loc[idx, s][::-1] == old.loc[idx, s]:
                    tmp = beagle_hla.loc[i, f'HLA-{gene} 1']
                    beagle_hla.loc[i, f'HLA-{gene} 1'] = beagle_hla.loc[i, f'HLA-{gene} 2']
                    beagle_hla.loc[i, f'HLA-{gene} 2'] = tmp
            else:
                beagle_hla.loc[i, f'HLA-{gene} 1'] = 'N/A'

        result = calculate_phasing_concordance(beagle_hla, our_hla, gene)
        result.to_csv(output.concordance_df, sep = '\t', header = True, index = False)

rule aggregate_phasing_concordance:
    input:
        concordance_df = expand("results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}-{gene}.tsv", gene = HLA_GENES, allow_missing = True)
    output:
        concordance_df = "results/phasing/oneKG_{vcf_version}-phasing-concordance-{filter}.tsv"
    resources: mem = '30G'
    threads: 1
    run: 
        df = pd.concat([pd.read_csv(f"results/phasing/oneKG_{wildcards.vcf_version}-phasing-concordance-{wildcards.filter}-{gene}.tsv", sep = '\t') for gene in HLA_GENES])

        df.to_csv(output.concordance_df, sep = '\t', header = True, index = False)

rule prepare_hla_reference_panel_method2:
    input:
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        reads1 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv"),
        reads2 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv"),
        mate_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv"),
        pair_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv")
    resources:
        mem = '120G'
    threads: 8
    params:
        outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}",
        hla_gene_information = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv",
        script = "/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/hla_align.py"
    shell: """
        python {params.script} {wildcards.hla_gene} {input.bam} {params.outdir}
    """

rule hla_imputation_method:
    input:
        bamlist = "results/hla/imputation/bamlists_fv/bamlist{num}.txt",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{hla_gene}fullallelesfilledin.RData",
        prepared_db = '/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv',
        reads1 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv", id = samples_fv, allow_missing = True),
        reads2 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv", id = samples_fv, allow_missing = True),
        mate_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv", id = samples_fv, allow_missing = True),
        pair_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv", id = samples_fv, allow_missing = True)
    output:
        imputed = "results/hla/imputation/QUILT_HLA_result_method/genes{num}/{hla_gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '60G'
    threads: 4
    params:
        quilt_sus_hla = tools['quilt_sus_hla'],
        fa_dict = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus1"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/genes{wildcards.num}/{wildcards.hla_gene}/

        {params.quilt_sus_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/genes{wildcards.num}/{wildcards.hla_gene}/" \
        --bamlist={input.bamlist} \
        --region={wildcards.hla_gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.hla_gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """

# oneKG_html = expand("results/phasing/html/oneKG-{filter}-{gene}.html", gene = HLA_GENES, filter = filters), 
# GAMCC_html = expand("results/phasing/html/GAMCC-{gene}.html", gene = HLA_GENES), 

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

rule prepare_hla_reference_panel_method2:
    input:
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        reads1 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv"),
        reads2 = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv"),
        mate_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv"),
        pair_matrix = temp("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv")
    resources:
        mem = '120G'
    threads: 8
    params:
        outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}",
        hla_gene_information = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv",
        script = "/well/band/users/rbx225/software/QUILT_sus/QUILT/Python/hla_align.py"
    shell: """
        python {params.script} {wildcards.hla_gene} {input.bam} {params.outdir}
    """

rule hla_imputation_prep_all:
    input:
        reads1 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv", hla_gene = HLA_GENES, id = samples_fv),
        reads2 = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv", hla_gene = HLA_GENES, id = samples_fv),
        mate_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv", hla_gene = HLA_GENES, id = samples_fv),
        pair_matrix = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv", hla_gene = HLA_GENES, id = samples_fv)

rule hla_imputation_prepare_align_per_allele:
    input:
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        read1_ary = temp("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary1.npy"),
        read2_ary = temp("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary2.npy")
    resources:
        mem = '40G'
    threads: 4
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv"
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        db = pd.read_csv(input.db_file, sep = ' ')
        gene = wildcards.hla_gene
        reads_apart_max = 10000
        bam = input.bam

        reads1 = get_chr6_reads(gene, bam, hla_gene_information, reads_apart_max)
        reads2 = get_hla_reads(gene, bam, reads_apart_max)

        if reads1.empty:
            reads1 = reads2.iloc[:2, :] if not reads2.empty else pd.DataFrame()
        elif reads2.empty:
            reads2 = reads1.iloc[:2, :]
        else:
            pass

        reads1['rev_seq'] = reads1['sequence'].apply(reverse_complement)
        reads1['rev_bq'] = reads1['base_quality'].apply(lambda bq: bq[::-1])

        j = wildcards.idx
        a = db.columns[j]
        res = per_allele(j, a, db, reads1)
        scores_ary = np.maximum(res[1], res[2])

        np.set_printoptions(precision=6, suppress=True)
        np.save(output.read1_ary, scores_ary)

        reads2['rev_seq'] = reads2['sequence'].apply(reverse_complement)
        reads2['rev_bq'] = reads2['base_quality'].apply(lambda bq: bq[::-1])

        res = per_allele(j, a, db, reads2)
        scores_ary = np.maximum(res[1], res[2])
        
        np.set_printoptions(precision=6, suppress=True)
        np.save(output.read2_ary, scores_ary)

rule hla_imputation_prepare_per_sample:
    input:
        read1_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary1.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        read2_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary2.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        bam = "data/bams/{id}.bam",
        db_file = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3570_aligners/{hla_gene}.ssv"
    output:
        reads1 = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads1.csv",
        reads2 = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/reads2.csv",
        mate_matrix = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/mate_likelihood_matrix.ssv",
        pair_matrix = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/alignment_likelihoods/{id}-{hla_gene}/pair_likelihood_matrix.ssv"
    resources:
        mem = '60G'
    threads: 4
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv"
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        db = pd.read_csv(input.db_file, sep = ' ')
        gene = wildcards.hla_gene
        reads_apart_max = 10000
        bam = input.bam

        reads1 = get_chr6_reads(gene, bam, hla_gene_information, reads_apart_max)
        reads2 = get_hla_reads(gene, bam, reads_apart_max)

        scores_mat1 = [np.hstack([np.load(f'results/hla/imputation/ref_panel/QUILT_bam_aligner/{wildcards.hla_gene}/{wildcards.id}/{wildcards.id}.{idx}.ary1.npy')[:,np.newaxis]]) for idx in range(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)])]
        likelihood_mat1 = np.exp(scores_mat1)/np.sum(np.exp(scores_mat1), axis = 1, keepdims = True)
        likemat1 = np.log(likelihood_mat1)

        scores_mat2 = [np.hstack([np.load(f'results/hla/imputation/ref_panel/QUILT_bam_aligner/{wildcards.hla_gene}/{wildcards.id}/{wildcards.id}.{idx}.ary2.npy')[:,np.newaxis]]) for idx in range(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)])]
        likelihood_mat2 = np.exp(scores_mat2)/np.sum(np.exp(scores_mat2), axis = 1, keepdims = True)
        likemat2 = np.log(likelihood_mat2)

        rl = reads1['sequence'].str.len().mode().values[0]
        n_mismatches = 5
        assumed_bq = 0.001
        min_valid_prob = np.log(math.comb(rl, n_mismatches)) + n_mismatches*np.log(assumed_bq) + (rl - n_mismatches)*np.log(1 - assumed_bq)

        valid_indices1 = np.any(likemat1 >= min_valid_prob, axis=1)
        valid_indices2 = np.any(likemat2 >= min_valid_prob, axis=1)
        likemat1, reads1 = likemat1[valid_indices1], reads1[valid_indices1]
        likemat2, reads2 = likemat2[valid_indices2], reads2[valid_indices2]
        
        likemat_all = np.vstack((likemat1, likemat2))

        id1, id2 = reads1.iloc[:, 0].to_numpy(), reads2.iloc[:, 0].to_numpy()
        
        readind = (reads1.iloc[:, 1].astype(int) // 64) % 4
        readind2 = (reads2.iloc[:, 1].astype(int) // 64) % 4
        mate_indicator = np.concatenate((readind, readind2))

        ids_all = np.concatenate((id1, id2))
        unique_ids = np.unique(ids_all)
        likemat_mate = np.zeros((len(unique_ids), likemat_all.shape[1]))

        for i, uid in enumerate(unique_ids):
            t1 = likemat_all[ids_all == uid, :]
            t2 = mate_indicator[ids_all == uid]
            if len(t2) > 0:
                likemat_mate[i, :] = np.sum(t1[t2 > 0], axis=0)

        valid_mask = likemat_mate.max(axis=1) >= min_valid_prob
        likemat_mate = likemat_mate[valid_mask]
        likemat_norm = 0.5 * np.exp(likemat_mate - likemat_mate.max(axis=1, keepdims=True)) + 1e-100

        likemat_paired = likemat_norm.T @ likemat_norm
        likemat_paired = pd.DataFrame(likemat_paired, index=db.columns, columns=db.columns)
        likemat_mate = pd.DataFrame(likemat_mate, index = unique_ids[valid_mask], columns=db.columns)

        reads1.to_csv(output.reads1, header = False, index = False)
        reads2.to_csv(output.reads2, header = False, index = False)
        pd.set_option('display.float_format', '{:.6e}'.format)
        likemat_mate.to_csv(output.mate_matrix, index=True, header=True, sep = ' ')
        likemat_paired.to_csv(output.pair_matrix, index=True, header=True, sep = ' ')

rule hla_la_calling:
    input:
        bam = "data/bams/{id}.bam",
        bai = "data/bams/{id}.bam.bai"
    output:
        called = "results/hla/call/{id}/hla/R1_bestguess_G.txt"
    resources:
        mem = '60G'
    threads: 4
    shell: """
        mkdir -p results/hla/call/{wildcards.id}/
        module load Java/17

        HLA-LA.pl \
        --BAM {input.bam} \
        --graph PRG_MHC_GRCh38_withIMGT \
        --workingDir /well/band/users/rbx225/GAMCC/results/hla/call/ \
        --sampleID {wildcards.id}
    """
called = expand("results/hla/call/{id}/hla/R1_bestguess_G.txt", id = samples_fv),

rule hla_alignment_new:
    input:
        read1_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary1.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        read2_ary = lambda wildcards: expand("results/hla/imputation/ref_panel/QUILT_bam_aligner/{hla_gene}/{id}/{id}.{idx}.ary2.npy", idx = np.arange(v3570_db_alleles[hla_genes.index(wildcards.hla_gene)]), allow_missing = True),
        bam = "data/bams/{id}.bam",
        db_files = expand("/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/{gene}.ssv", gene = HLA_GENES_ALL_EXPANDED)
    output:
        reads1 = temp("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/reads1.csv"),
        reads2 = temp("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/reads2.csv"),
        mate_matrix = temp("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/mate_likelihood_matrix.ssv"),
        pair_matrix = temp("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/pair_likelihood_matrix.ssv"),
        read_topresult = "results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.onlyreads.topresult.txt"
    resources:
        mem = '60G'
    threads: 4
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv",
        db_dir = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_oneKG_only/",
        tmp_outdir = "results/hla/imputation/QUILT_HLA_result_method/",
        reads_df_outdir = "results/hla/imputation/QUILT_HLA_result_method/"
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        reads_apart_max = 1000
        reads_extend_max = 1000
        
        reads1 = get_chr6_reads(wildcards.gene, input.bam, hla_gene_information, 
                        reads_apart_max = reads_apart_max, 
                        reads_extend_max = reads_extend_max)
    
        rl = reads1['sequence'].str.len().mode().values[0]

        reads2 = get_hla_reads(wildcards.gene, input.bam, reads_extend_max = reads_extend_max)

        if reads1.empty:
            reads1 = reads2.iloc[:2, :] if not reads2.empty else pd.DataFrame()
        elif reads2.empty:
            reads2 = reads1.iloc[:2, :]
        else:
            pass
        
        ncores = 2*(len(os.sched_getaffinity(0))) - 1

        db_dict = {}
        for g in HLA_GENES_ALL:
            db = pd.read_csv(f'{params.db_dir}{g}.ssv', sep = ' ')
            for c in db.columns:
                refseq = ''.join(db[c].tolist()).replace('.', '').lstrip('*').rstrip('*')
                db_dict[c] = refseq

        columns = np.array(list(db_dict.keys()))
        likemat1 = multi_calculate_loglikelihood_per_allele(reads1, db_dict, ncores)
        min_valid_prob = n_mismatches*np.log(assumed_bq) + (rl - n_mismatches)*np.log(1 - assumed_bq)

        valid_indices1 = np.any(likemat1 >= min_valid_prob, axis=1)
        likemat1, reads1 = likemat1[valid_indices1], reads1[valid_indices1]
        
        likemat_all = likemat1
        id1 = reads1.iloc[:, 0].to_numpy()

        unique_ids = np.unique(id1)

        for g in HLA_GENES:
        likemat_mate = -600*np.ones((len(unique_ids), likemat_all.shape[1]))

        for i, uid in enumerate(unique_ids):
            tmp = likemat_all[id1 == uid, :]

            keep = False
            for j in range(tmp.shape[0]):
                best_indices = np.where(tmp[j,:] == tmp[j,:].max())[0]
                best_alleles = np.array(columns)[best_indices]
                aligned_genes = np.unique(np.array([s.split('*')[0] for s in best_alleles]))
                keep = keep or ((len(aligned_genes) == 1) and (aligned_genes[0] == gene))

            if keep:
                likemat_mate[i, :] = np.sum(tmp, axis=0)

        columns_to_keep = np.where(np.char.startswith(columns, gene + '*'))[0]
        scoring_df = pd.DataFrame({'ID': unique_ids, 
                            'Target': likemat_mate[:,columns_to_keep].max(axis = 1), 
                            'Others': likemat_mate[:,~np.where(np.char.startswith(columns, gene + '*'))[0]].max(axis = 1)
                        })
        scoring_df['diff'] = scoring_df['Target'] - scoring_df['Others']
        valid_mask = ((scoring_df['diff'] >= score_diff_in_alignment_genes) & (scoring_df['Target'] >= min_valid_prob)).tolist()

        if len(np.where(valid_mask)[0]) <= 2:
            valid_mask = (scoring_df['Target'] >= min_valid_prob).tolist()

        scoring_df = scoring_df[valid_mask]
        target_alleles = columns[columns_to_keep]
        likemat_mate = likemat_mate[valid_mask, :][:, columns_to_keep]
        unique_ids = unique_ids[valid_mask]
        reads1 = reads1[reads1['ID'].isin(unique_ids)]
        likemat_mate_df = pd.DataFrame(likemat_mate, index = unique_ids, columns=target_alleles)

        reads1_archive = reads1.copy()
        likemat_mate_df_archive = likemat_mate_df.copy()
        likemat_mate_archive = likemat_mate.copy()

        if likemat_mate.shape[0] > 2:
            # Filtering out reads aligned to more than one-field resolution
            df = likemat_mate_df
            best_alignments = df.max(axis=1)
            matching_cols = df.eq(best_alignments, axis=0).apply(lambda row: likemat_mate_df.columns[row].tolist(), axis=1)
            one_field = matching_cols.apply(lambda cols: list(set(col.split(':')[0] for col in cols)))
            counts = one_field.apply(len)
            rows_to_remove = counts[counts >= n_unique_onefield].index
            rows_to_keep = df.index.difference(rows_to_remove)
            likemat_mate_df = likemat_mate_df.loc[rows_to_keep]
            reads1 = reads1[reads1['ID'].isin(rows_to_keep)]
            likemat_mate = likemat_mate_df.values
            unique_ids = reads1['ID'].tolist()

            if likemat_mate.shape[0] < 2:
                reads1 = reads1_archive.copy()
                likemat_mate_df = likemat_mate_df_archive.copy()
                likemat_mate = likemat_mate_archive.copy()
                unique_ids = reads1['ID'].unique().tolist()

        if likemat_mate.shape[0] > 2:
            likemat_mate_normalised = normalize(np.exp(likemat_mate - likemat_mate.max(axis = 1, keepdims = True)), axis=1, norm='l2')
            gmm = GaussianMixture(n_components=2, tol = 1e-10, max_iter = 1000, n_init = 10)
            gmm.fit(likemat_mate_normalised)
            labels = gmm.predict(likemat_mate_normalised)

            group1 = np.where(labels == 0)[0]
            group2 = np.where(labels == 1)[0]

            l1 = likemat_mate[group1, :].sum(axis=0)
            l2 = likemat_mate[group2, :].sum(axis=0)

            likemat_mate_phased = np.vstack([l1, l2])
        else:
            likemat_mate_phased = likemat_mate

        likemat_norm=likemat_mate_phased-likemat_mate_phased.max(axis = 1, keepdims = True)
        likemat_norm=0.5*np.exp(likemat_norm)+1e-100

        likemat_pair=pd.DataFrame(0, index=target_alleles, columns=target_alleles)
        qq=likemat_pair*0
        for i in range(likemat_mate_phased.shape[0]):
            qq=qq*0
            m1=qq+likemat_norm[i,:]
            m1=m1+m1.T
            likemat_pair=likemat_pair+np.log(m1)

        likemat_paired_df = pd.DataFrame(likemat_pair, index=target_alleles, columns=target_alleles)

        reads1.to_csv(output.reads1, header = False, index = False)
        reads2.to_csv(output.reads2, header = False, index = False)
        likemat_mate_df.to_csv(output.mate_matrix,  index=True, header=True, sep = ' ')
        likemat_paired_df.to_csv(output.pair_matrix,  index=True, header=True, sep = ' ')
        
        logllmax = likemat_paired_df.max().max()
        row_index, col_index = np.where(likemat_paired_df == logllmax)
        a1, a2 = likemat_paired_df.index[row_index[0]], likemat_paired_df.columns[col_index[0]]

        if a1.count(':') > 1:
            a1 = ':'.join(a1.split(':')[:2])
        if a2.count(':') > 1:
            a2 = ':'.join(a2.split(':')[:2])

        result = pd.DataFrame({'sample_number': 1,
                            'sample_name': wildcards.id,
                            'bestallele1': a1,
                            'bestallele2': a2,
                            'post_prob': logllmax,
                            'sums': 0}, index = ['0'])
        result.to_csv(output.read_topresult, index = False, header = True, sep = '\t')

rule hla_alignment_matrix:
    input:
        bam = "data/bams/{id}.bam",
        db_files = expand("/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_aligners/{gene}.ssv", gene = HLA_GENES_ALL_EXPANDED)
    output:
        matrix = "results/hla/imputation/WFA_alignments/v3390_oneKG/{id}/{gene}/AS_matrix.ssv"
    resources:
        mem = '60G'
    threads: 8
    params:
        hla_gene_information_file = "/well/band/users/rbx225/software/QUILT_sus/hla_ancillary_files/hla_gene_information.tsv",
        db_dir = "/well/band/users/rbx225/recyclable_files/hla_reference_files/v3390_oneKG_only/",
        tmp_outdir = "results/hla/imputation/QUILT_HLA_result_method/",
    run:
        hla_gene_information = pd.read_csv(params.hla_gene_information_file, sep = ' ')
        reads_apart_max = 1000
        reads_extend_max = 1000
        
        reads1 = get_chr6_reads(wildcards.gene, input.bam, hla_gene_information, 
                        reads_apart_max = reads_apart_max, 
                        reads_extend_max = reads_extend_max)
    
        rl = reads1['sequence'].str.len().mode().values[0]
        ncores = 7

        db_dict = {}
        for g in HLA_GENES_ALL:
            db = pd.read_csv(f'{params.db_dir}{g}.ssv', sep = ' ')
            for c in db.columns:
                refseq = ''.join(db[c].tolist()).replace('.', '').lstrip('*').rstrip('*')
                db_dict[c] = refseq

        columns = np.array(list(db_dict.keys()))
        likemat1 = multi_calculate_loglikelihood_per_allele(reads1, db_dict, ncores)
        
        alignment_scores_raw = pd.DataFrame(likemat1, index=reads1['ID'].tolist(), columns=columns)
        alignment_scores_raw.to_csv(output.matrix, float_format='%.6f', sep = ' ', header = True, index = True)
        
# ncores = 2*(len(os.sched_getaffinity(0))) - 1

rule all:
    input:
        filtered_vcf = expand("results/wip_vcfs/malariaGen_v1_b38/vanilla/topmed_sites/lc.chr{chr}.vcf.gz", chr = chromosome)

rule get_TOPMed_sites:
    input:
        vcf = "results/two-stage-imputation/vanilla/malariaGen_v1_b38_topmed/vcf/chr{chr}.dose.vcf.gz"
    output:
        sites = "data/topmed_variants/chr{chr}.tsv"
    resources:
        mem = '40G'
    threads: 4
    shell: """
        bcftools query -f '%CHROM\t%POS\n' {input.vcf} | sort -u > {output.sites}
    """

rule filter_topmed_sites:
    input:
        vcf = "results/wip_vcfs/malariaGen_v1_b38/vanilla/high_info_high_af/lc.chr{chr}.vcf.gz",
        sites = "data/topmed_variants/chr{chr}.tsv"
    output:
        filtered_vcf = "results/wip_vcfs/malariaGen_v1_b38/vanilla/topmed_sites/lc.chr{chr}.vcf.gz"
    resources:
        mem = '40G'
    threads: 4
    shell: """
        mkdir -p results/wip_vcfs/malariaGen_v1_b38/vanilla/topmed_sites/

        bcftools view -R {input.sites} -Oz -o {output.filtered_vcf} {input.vcf}
        tabix -f {output.filtered_vcf}
    """