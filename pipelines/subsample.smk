configfile: "pipelines/config.json"

clean_fastq = config['clean_fastq']
concatenate = config['concatenate']
subsample_depth = config['subsample_depth_1x']*config['subsample_depth']

rule ss_fastq:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz" if clean_fastq else "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz" if clean_fastq else "data/fastq/{id}_2.fastq.gz"
    output:
        ss_fastq1 = "data/subsampled_fastq/{id}_subsampled_1.fastq",
        ss_fastq2 = "data/subsampled_fastq/{id}_subsampled_2.fastq"
    resources:
        mem_mb = 30000
    params: 
        n = subsample_depth
    shell: """
        seqtk sample -s100 {input.fastq1} {params.n} > {output.ss_fastq1}
        seqtk sample -s100 {input.fastq2} {params.n} > {output.ss_fastq2}
    """

rule ss_alignment:
    input:
        ss_fastq1 = "data/subsampled_fastq/{id}_subsampled_1.fastq",
        ss_fastq2 = "data/subsampled_fastq/{id}_subsampled_2.fastq",
        reference = "data/references/concatenated/GRCh38_no_alt_Pf3D7_v3_phiX.fasta" if concatenate else config["ref38"]
    output:
        bam = temp("data/subsampled_bams/tmp/{id}_subsampled.bam")
    resources:
        mem_mb = 30000
    threads: 1
    shell: """
	    bwa mem -t {threads} {input.reference} {input.ss_fastq1} {input.ss_fastq2} | samtools view -b -o {output.bam}
    """

rule ss_fixmate:
	input:
		bam = rules.ss_alignment.output.bam
	output:
		fixmate = temp("data/subsampled_bams/tmp/{id}_subsampled.fixmate.bam")
	resources: mem_mb = 30000
	shell: """
		samtools fixmate -m {input.bam} {output.fixmate}
	"""

rule ss_sort:
	input:
		fixmate = rules.ss_fixmate.output.fixmate
	output:
		sorted = temp("data/subsampled_bams/tmp/{id}_subsampled.sorted.bam")
	resources: mem_mb = 30000
	shell: """
		samtools sort -o {output.sorted} {input.fixmate}
	"""

rule ss_index:
	input:
		sorted = rules.ss_sort.output.sorted
	output:
		bai = temp("data/subsampled_bams/tmp/{id}_subsampled.sorted.bam.bai")
	resources: mem_mb = 30000
	shell: """
		samtools index {input.sorted}
	"""

rule ss_rename_alignment_files:
	input:
		bam = rules.ss_index.input.sorted,
		bai = rules.ss_index.output.bai
	output:
		bam = "data/subsampled_bams/{id}_subsampled.bam",
		bai = "data/subsampled_bams/{id}_subsampled.bam.bai"
	shell: """
		mv {input.bam} {output.bam}
		mv {input.bai} {output.bai}
	"""
