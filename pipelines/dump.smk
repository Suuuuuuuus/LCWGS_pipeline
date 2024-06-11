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