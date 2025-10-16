configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"

subsample_depth = int(config['subsample_depth_1x']*config['subsample_depth'])
clean_fastq = config['clean_fastq']

if config['clean_fastq']:
    ruleorder: ss_fastq_alt > ss_fastq
else:
    ruleorder: ss_fastq > ss_fastq_alt

rule ss_fastq:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        ss_fastq1 = temp("data/subsampled_fastq/{id}_subsampled_1.fastq"),
        ss_fastq2 = temp("data/subsampled_fastq/{id}_subsampled_2.fastq")
    resources:
        mem = '30G'
    params:
        n = subsample_depth
    shell: """
        zcat {input.fastq1} | seqtk sample -s100 - {params.n} > {output.ss_fastq1}
        zcat {input.fastq2} | seqtk sample -s100 - {params.n} > {output.ss_fastq2}
    """

rule ss_fastq_alt:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz"
    output:
        ss_fastq1 = temp("data/subsampled_fastq/{id}_subsampled_1.fastq"),
        ss_fastq2 = temp("data/subsampled_fastq/{id}_subsampled_2.fastq")
    resources:
        mem = '30G'
    params:
        n = subsample_depth
    shell: """
        zcat {input.fastq1} | seqtk sample -s100 - {params.n} > {output.ss_fastq1}
        zcat {input.fastq2} | seqtk sample -s100 - {params.n} > {output.ss_fastq2}
    """

rule ss_alignment:
    input:
        ss_fastq1 = temp("data/subsampled_fastq/{id}_subsampled_1.fastq"),
        ss_fastq2 = temp("data/subsampled_fastq/{id}_subsampled_2.fastq"),
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta"
    output:
        bam = "data/subsampled_bams/{id}_subsampled.bam",
        bai = "data/subsampled_bams/{id}_subsampled.bam.bai",
        tmp1 = temp("data/subsampled_bams/{id}_subsampled.tmp1.bam"),
        metric = temp("data/subsampled_bams/{id}_subsampled.metrics.txt")
    resources:
        mem = '30G'
    params: 
        sample = "{id}",
        picard = tools["picard_plus"]
    threads: 6
    shell: """
        bwa mem -t {threads} {input.reference} {input.ss_fastq1} {input.ss_fastq2} | samtools view -b -o {output.tmp1}
        
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